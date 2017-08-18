use std::cmp::Ordering;
use std::collections::BinaryHeap;
use num_traits::Float;
use types::{Point, LineString, Polygon, MultiLineString, MultiPolygon};

use spade::SpadeFloat;
use spade::primitives::SimpleEdge;
use spade::BoundingRect;
use spade::rtree::RTree;

// Store triangle information
// current is the candidate point for removal
#[derive(Debug)]
struct VScore<T>
where
    T: Float,
{
    left: usize,
    current: usize,
    right: usize,
    area: T,
}

// These impls give us a min-heap
impl<T> Ord for VScore<T>
where
    T: Float,
{
    fn cmp(&self, other: &VScore<T>) -> Ordering {
        other.area.partial_cmp(&self.area).unwrap()
    }
}

impl<T> PartialOrd for VScore<T>
where
    T: Float,
{
    fn partial_cmp(&self, other: &VScore<T>) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl<T> Eq for VScore<T>
where
    T: Float,
{
}

impl<T> PartialEq for VScore<T>
where
    T: Float,
{
    fn eq(&self, other: &VScore<T>) -> bool
    where
        T: Float,
    {
        self.area == other.area
    }
}

// Simplify a line using the [Visvalingam-Whyatt](http://www.tandfonline.com/doi/abs/10.1179/000870493786962263) algorithm
//
// epsilon is the minimum triangle area
// The paper states that:
// If [the new triangle's] calculated area is less than that of the last point to be
// eliminated, use the latter's area instead.
// (This ensures that the current point cannot be eliminated
// without eliminating previously eliminated points)
// (Visvalingam and Whyatt 2013, p47)
// However, this does *not* apply if you're using a user-defined epsilon;
// It's OK to remove triangles with areas below the epsilon,
// then recalculate the new triangle area and push it onto the heap
// based on Huon Wilson's original implementation:
// https://github.com/huonw/isrustfastyet/blob/25e7a68ff26673a8556b170d3c9af52e1c818288/mem/line_simplify.rs
fn visvalingam<T>(orig: &[Point<T>], epsilon: &T) -> Vec<Point<T>>
where
    T: Float,
{
    // No need to continue without at least three points
    if orig.len() < 3 || orig.is_empty() {
        return orig.to_vec();
    }

    let max = orig.len();

    // Adjacent retained points. Simulating the points in a
    // linked list with indices into `orig`. Big number (larger than or equal to
    // `max`) means no next element, and (0, 0) means deleted element.
    let mut adjacent: Vec<(_)> = (0..orig.len())
        .map(|i| if i == 0 {
            (-1_i32, 1_i32)
        } else {
            ((i - 1) as i32, (i + 1) as i32)
        })
        .collect();

    // Store all the triangles in a minimum priority queue, based on their area.
    // Invalid triangles are *not* removed if / when points
    // are removed; they're handled by skipping them as
    // necessary in the main loop by checking the corresponding entry in
    // adjacent for (0, 0) values
    let mut pq = BinaryHeap::new();
    // Compute the initial triangles, i.e. take all consecutive groups
    // of 3 points and form triangles from them
    for (i, win) in orig.windows(3).enumerate() {
        pq.push(VScore {
            area: area(win.first().unwrap(), &win[1], win.last().unwrap()),
            current: i + 1,
            left: i,
            right: i + 2,
        });
    }
    // While there are still points for which the associated triangle
    // has an area below the epsilon
    loop {
        let smallest = match pq.pop() {
            // We've exhausted all the possible triangles, so leave the main loop
            None => break,
            // This triangle's area is above epsilon, so skip it
            Some(ref x) if x.area > *epsilon => continue,
            //  This triangle's area is below epsilon: eliminate the associated point
            Some(s) => s,
        };
        let (left, right) = adjacent[smallest.current];
        // A point in this triangle has been removed since this VScore
        // was created, so skip it
        if left as i32 != smallest.left as i32 || right as i32 != smallest.right as i32 {
            continue;
        }
        // We've got a valid triangle, and its area is smaller than epsilon, so
        // remove it from the simulated "linked list"
        let (ll, _) = adjacent[left as usize];
        let (_, rr) = adjacent[right as usize];
        adjacent[left as usize] = (ll, right);
        adjacent[right as usize] = (left, rr);
        adjacent[smallest.current as usize] = (0, 0);

        // Now recompute the triangle area, using left and right adjacent points
        let choices = [(ll, left, right), (left, right, rr)];
        for &(ai, current_point, bi) in &choices {
            if ai as usize >= max || bi as usize >= max {
                // Out of bounds, i.e. we're on one edge
                continue;
            }
            let new_left = Point::new(orig[ai as usize].x(), orig[ai as usize].y());
            let new_current = Point::new(
                orig[current_point as usize].x(),
                orig[current_point as usize].y(),
            );
            let new_right = Point::new(orig[bi as usize].x(), orig[bi as usize].y());
            pq.push(VScore {
                area: area(&new_left, &new_current, &new_right),
                current: current_point as usize,
                left: ai as usize,
                right: bi as usize,
            });
        }
    }
    // Filter out the points that have been deleted, returning remaining points
    orig.iter()
        .zip(adjacent.iter())
        .filter_map(|(tup, adj)| if *adj != (0, 0) { Some(*tup) } else { None })
        .collect::<Vec<Point<T>>>()
}

// Visvalingam with self-intersection detection to preserve topologies
// this is a port of the technique at https://www.jasondavies.com/simplify/
fn visvalingam_preserve<T>(orig: &[Point<T>], epsilon: &T) -> Vec<Point<T>>
where
    T: Float + SpadeFloat,
{
    let mut internal_epsilon = *epsilon;
    // No need to continue without at least three points
    if orig.len() < 3 || orig.is_empty() {
        return orig.to_vec();
    }

    let max = orig.len();

    // Adjacent retained points. Simulating the points in a
    // linked list with indices into `orig`. Big number (larger than or equal to
    // `max`) means no next element, and (0, 0) means deleted element.
    let mut adjacent: Vec<(_)> = (0..orig.len())
        .map(|i| if i == 0 {
            (-1_i32, 1_i32)
        } else {
            ((i - 1) as i32, (i + 1) as i32)
        })
        .collect();
    let mut intersections = vec![];

    // Store all the triangles in a minimum priority queue, based on their area.
    // Invalid triangles are *not* removed if / when points
    // are removed; they're handled by skipping them as
    // necessary in the main loop by checking the corresponding entry in
    // adjacent for (0, 0) values
    let mut pq = BinaryHeap::new();
    let mut tree: RTree<SimpleEdge<_>> = RTree::new();
    // Compute the initial triangles, i.e. take all consecutive groups
    // of 3 points and form triangles from them
    for (i, win) in orig.windows(3).enumerate() {
        let v = VScore {
            area: area(&win[0], &win[1], &win[2]),
            current: i + 1,
            left: i,
            right: i + 2,
        };
        pq.push(v);
        // populate R* tree with line segments
        tree.insert(SimpleEdge::new(win[0], win[1]));
        tree.insert(SimpleEdge::new(win[1], win[2]));
    }
    // While there are still points for which the associated triangle
    // has an area below the epsilon
    loop {
        let mut smallest = match pq.pop() {
            // We've exhausted all the possible triangles, so leave the main loop
            None => break,
            // This triangle's area is above epsilon, so skip it
            Some(ref x) if x.area > internal_epsilon => continue,
            //  This triangle's area is below epsilon: eliminate the associated point
            Some(s) => s,
        };
        let (left, right) = adjacent[smallest.current];
        // A point in this triangle has been removed since this VScore
        // was created, so skip it
        if left as i32 != smallest.left as i32 || right as i32 != smallest.right as i32 {
            continue;
        }
        // if removal of this point causes an intersection, save it for re-processing
        if tree_intersect(&tree, &mut smallest, orig) {
            // TODO: ensure that we're doing this at the correct point in the loop
            // decrease area of next-largest triangle in heap to epsilon
            // this means that we remove smallest and the next-largest in order
            let mut to_alter = pq.peek_mut().unwrap();
            to_alter.area = internal_epsilon;
            &intersections.push(smallest);
            continue;
        }
        while !intersections.is_empty() {
            pq.push(intersections.pop().unwrap());
        }
        // We've got a valid triangle, and its area is smaller than epsilon, so
        // remove it from the simulated "linked list"
        adjacent[smallest.current as usize] = (0, 0);
        // Now recompute the triangle area, using left and right adjacent points
        let (ll, _) = adjacent[left as usize];
        let (_, rr) = adjacent[right as usize];
        adjacent[left as usize] = (ll, right);
        adjacent[right as usize] = (left, rr);
        let choices = [(ll, left, right), (left, right, rr)];
        for &(ai, current_point, bi) in &choices {
            if ai as usize >= max || bi as usize >= max {
                // Out of bounds, i.e. we're on one edge
                continue;
            }
            let new_left = Point::new(orig[ai as usize].x(), orig[ai as usize].y());
            let new_current = Point::new(
                orig[current_point as usize].x(),
                orig[current_point as usize].y(),
            );
            let new_right = Point::new(orig[bi as usize].x(), orig[bi as usize].y());
            let new_triangle = VScore {
                area: area(&new_left, &new_current, &new_right),
                current: current_point as usize,
                left: ai as usize,
                right: bi as usize,
            };
            // we have to call this twice because only one segment is returned at a time
            // this should be OK because a point can only share at most two segments
            tree.lookup_and_remove(&orig[smallest.right]);
            tree.lookup_and_remove(&orig[smallest.left]);
            // add re-computed line segments to the tree
            tree.insert(SimpleEdge::new(orig[ai as usize], orig[current_point as usize]));
            tree.insert(SimpleEdge::new(orig[current_point as usize], orig[bi as usize]));
            // push re-computed triangle onto heap
            pq.push(new_triangle);
        }
    }
    // Filter out the points that have been deleted, returning remaining points
    orig.iter()
        .zip(adjacent.iter())
        .filter_map(|(tup, adj)| if *adj != (0, 0) { Some(*tup) } else { None })
        .collect::<Vec<Point<T>>>()
}

// is p1 -> p2 -> p3 wound counterclockwise?
fn ccw<T>(p1: &Point<T>, p2: &Point<T>, p3: &Point<T>) -> bool
where
    T: Float,
{
    (p3.y() - p1.y()) * (p2.x() - p1.x()) > (p2.y() - p1.y()) * (p3.x() - p1.x())
}

// checks whether line segments with p1-p4 as their start and endpoints touch or cross
fn cartesian_intersect<T>(p1: &Point<T>, p2: &Point<T>, p3: &Point<T>, p4: &Point<T>) -> bool
where
    T: Float,
{
    (ccw(p1, p3, p4) ^ ccw(p2, p3, p4)) & (ccw(p1, p2, p3) ^ ccw(p1, p2, p4))
}

// check whether a triangle's edges intersect with any other edges of the LineString
fn tree_intersect<T>(tree: &RTree<SimpleEdge<Point<T>>>, triangle: &VScore<T>, orig: &[Point<T>]) -> bool
where
    T: Float + SpadeFloat,
{
    let point_a = orig[triangle.left];
    let point_c = orig[triangle.right];
    let candidates = tree.lookup_in_rectangle(&BoundingRect::from_corners(&point_a, &point_c));
    candidates
        .iter()
        .map(|c| {
            // triangle start point, end point
            let ca = c.from;
            let cb = c.to;
            if ca != point_a && ca != point_c && cb != point_a && cb != point_c && cartesian_intersect(&ca, &cb, &point_a, &point_c) {
                true
            } else {
                false
            }
        })
        .any(|elem| elem == true)
}

// Area of a triangle given three vertices
fn area<T>(p1: &Point<T>, p2: &Point<T>, p3: &Point<T>) -> T
where
    T: Float,
{
    ((p1.x() - p3.x()) * (p2.y() - p3.y()) - (p2.x() - p3.x()) * (p1.y() - p3.y())).abs() / (T::one() + T::one())
}

/// Simplifies a geometry.
///
/// Polygons are simplified by running the algorithm on all their constituent rings.  This may
/// result in invalid Polygons, and has no guarantee of perserving topology.  Multi* objects are
/// simplified by simplifying all their constituent geometries individually.
pub trait SimplifyVW<T, Epsilon = T> {
    /// Returns the simplified representation of a geometry, using the [Visvalingam-Whyatt](http://www.tandfonline.com/doi/abs/10.1179/000870493786962263) algorithm
    ///
    /// See [here](https://bost.ocks.org/mike/simplify/) for a graphical explanation
    ///
    /// ```
    /// use geo::{Point, LineString};
    /// use geo::algorithm::simplifyvw::{SimplifyVW};
    ///
    /// let mut vec = Vec::new();
    /// vec.push(Point::new(5.0, 2.0));
    /// vec.push(Point::new(3.0, 8.0));
    /// vec.push(Point::new(6.0, 20.0));
    /// vec.push(Point::new(7.0, 25.0));
    /// vec.push(Point::new(10.0, 10.0));
    /// let linestring = LineString(vec);
    /// let mut compare = Vec::new();
    /// compare.push(Point::new(5.0, 2.0));
    /// compare.push(Point::new(7.0, 25.0));
    /// compare.push(Point::new(10.0, 10.0));
    /// let ls_compare = LineString(compare);
    /// let simplified = linestring.simplifyvw(&30.0);
    /// assert_eq!(simplified, ls_compare)
    /// ```
    fn simplifyvw(&self, epsilon: &T) -> Self
    where
        T: Float;
}

/// Simplifies a geometry, preserving its topology by removing self-intersections
pub trait SimplifyVWPreserve<T, Epsilon = T> {
    /// Returns the simplified representation of a geometry, using the [Visvalingam-Whyatt](http://www.tandfonline.com/doi/abs/10.1179/000870493786962263) algorithm
    ///
    /// See [here](https://www.jasondavies.com/simplify/) for a graphical explanation
    ///
    /// The topology-preserving algorithm uses an [R* tree](../../../spade/rtree/struct.RTree.html) to efficiently find candidate line segments
    /// which are tested for intersection with a given triangle. If intersections are found,
    /// the triangle is stored for later processing, and the next-largest area in the minimum
    /// priority queue is chosen as the new epsilon.
    ///
    /// In the example below, `(135.0, 68.0)` would be retained by the standard algorithm,
    /// thus causing a self-intersection when the given epsilon is used. By increasing the epsilon
    /// to the area associated with that point, we ensure that it is removed, thus removing the self-intersection.
    ///
    /// ```
    /// use geo::{Point, LineString};
    /// use geo::algorithm::simplifyvw::{SimplifyVWPreserve};
    ///
    /// let mut vec = Vec::new();
    /// vec.push(Point::new(10., 60.));
    /// vec.push(Point::new(135., 68.));
    /// vec.push(Point::new(94., 48.));
    /// vec.push(Point::new(126., 31.));
    /// vec.push(Point::new(280., 19.));
    /// vec.push(Point::new(117., 48.));
    /// vec.push(Point::new(300., 40.));
    /// vec.push(Point::new(301., 10.));
    /// let linestring = LineString(vec);
    /// let mut compare = Vec::new();
    /// compare.push(Point::new(10., 60.));
    /// compare.push(Point::new(126., 31.));
    /// compare.push(Point::new(280., 19.));
    /// compare.push(Point::new(117., 48.));
    /// compare.push(Point::new(300., 40.));
    /// compare.push(Point::new(301., 10.));
    /// let ls_compare = LineString(compare);
    /// let simplified = linestring.simplifyvw_preserve(&668.6);
    /// assert_eq!(simplified, ls_compare)
    /// ```
    fn simplifyvw_preserve(&self, epsilon: &T) -> Self
    where
        T: Float + SpadeFloat;
}

impl<T> SimplifyVWPreserve<T> for LineString<T>
where
    T: Float + SpadeFloat,
{
    fn simplifyvw_preserve(&self, epsilon: &T) -> LineString<T> {
        LineString(visvalingam_preserve(&self.0, epsilon))
    }
}

impl<T> SimplifyVW<T> for LineString<T>
where
    T: Float,
{
    fn simplifyvw(&self, epsilon: &T) -> LineString<T> {
        LineString(visvalingam(&self.0, epsilon))
    }
}

impl<T> SimplifyVW<T> for MultiLineString<T>
where
    T: Float,
{
    fn simplifyvw(&self, epsilon: &T) -> MultiLineString<T> {
        MultiLineString(self.0.iter().map(|l| l.simplifyvw(epsilon)).collect())
    }
}

impl<T> SimplifyVW<T> for Polygon<T>
where
    T: Float,
{
    fn simplifyvw(&self, epsilon: &T) -> Polygon<T> {
        Polygon::new(
            self.exterior.simplifyvw(epsilon),
            self.interiors
                .iter()
                .map(|l| l.simplifyvw(epsilon))
                .collect(),
        )
    }
}


impl<T> SimplifyVW<T> for MultiPolygon<T>
where
    T: Float,
{
    fn simplifyvw(&self, epsilon: &T) -> MultiPolygon<T> {
        MultiPolygon(self.0.iter().map(|p| p.simplifyvw(epsilon)).collect())
    }
}

#[cfg(test)]
mod test {
    use types::{Point, LineString, Polygon, MultiLineString, MultiPolygon};
    use super::{visvalingam, visvalingam_preserve, SimplifyVW};
    use std::fs::File;
    use std::io::{Write, BufWriter};
    use num_traits::ToPrimitive;
    extern crate serde_json;
    use self::serde_json::{Map, to_value};
    extern crate geojson;
    use self::geojson::{Feature, GeoJson, Geometry, Value};

    #[test]
    fn visvalingam_test() {
        // this is the PostGIS example
        let points = vec![
            (5.0, 2.0),
            (3.0, 8.0),
            (6.0, 20.0),
            (7.0, 25.0),
            (10.0, 10.0),
        ];
        let points_ls: Vec<_> = points.iter().map(|e| Point::new(e.0, e.1)).collect();

        let correct = vec![(5.0, 2.0), (7.0, 25.0), (10.0, 10.0)];
        let correct_ls: Vec<_> = correct.iter().map(|e| Point::new(e.0, e.1)).collect();

        let simplified = visvalingam(&points_ls, &30.);
        assert_eq!(simplified, correct_ls);
    }
    #[test]
    fn quux() {
        // this LineString will have a self-intersection if the point with the
        // smallest associated area is removed
        // the associated triangle is (1, 2, 3), and has an area of 668.5
        // the new triangle (0, 1, 3) self-intersects with triangle (3, 4, 5)
        // by detecting the intersection, peeking at the next-largest value in the queue,
        // and increasing the epsilon, we ensure that
        // (0, 1, 2) is removed before re-processing (1, 2, 3)
        // our final LineString is (0, 3, 4, 5, 6, 7)
        let points = vec![
            (10., 60.),
            (135., 68.),
            (94., 48.),
            (126., 31.),
            (280., 19.),
            (117., 48.),
            (300., 40.),
            (301., 10.),
        ];
        let points_ls: Vec<_> = points.iter().map(|e| Point::new(e.0, e.1)).collect();
        let simplified = visvalingam_preserve(&points_ls, &668.6);
        // this is the correct, non-intersecting LineString
        let correct = vec![
            (10., 60.),
            (126., 31.),
            (280., 19.),
            (117., 48.),
            (300., 40.),
            (301., 10.),
        ];
        let correct_ls: Vec<_> = correct.iter().map(|e| Point::new(e.0, e.1)).collect();
        assert_eq!(simplified, correct_ls);
    }
    #[test]
    fn grim() {
        // simplify a longer LineString, hopefully eliminating self-intersections
        let points = include!("test_fixtures/norway_main.rs");
        let points_ls: Vec<_> = points.iter().map(|e| Point::new(e[0], e[1])).collect();
        let simplified = visvalingam_preserve(&points_ls, &0.0005);
        // dump simplified geometry to GeoJSON so it can be checked for validity elsewhere
        // let simplified_vec: Vec<Vec<f64>> = simplified
        //     .iter()
        //     .map(|point| {
        //         vec![point.x().to_f64().unwrap(), point.y().to_f64().unwrap()]
        //     })
        //     .collect();
        // let mut properties = Map::new();
        // properties.insert(
        //     String::from("name"),
        //     to_value("Norway_VW_Decrease").unwrap(),
        // );
        // let geometry = Geometry::new(Value::Polygon(vec![simplified_vec]));
        // let geojson = GeoJson::Feature(Feature {
        //     bbox: None,
        //     geometry: Some(geometry),
        //     id: None,
        //     properties: Some(properties),
        //     foreign_members: None,
        // });
        // let geojson_string = geojson.to_string();
        // let f = File::create("src/algorithm/test_fixtures/norway_reduced.geojson").unwrap();
        // let mut bw = BufWriter::new(f);
        // bw.write_all(geojson_string.as_bytes()).unwrap();
        assert_eq!(simplified.len(), 3267);
    }
    #[test]
    fn visvalingam_test_long() {
        // simplify a longer LineString
        let points = include!("test_fixtures/vw_orig.rs");
        let points_ls: Vec<_> = points.iter().map(|e| Point::new(e[0], e[1])).collect();
        let correct = include!("test_fixtures/vw_simplified.rs");
        let correct_ls: Vec<_> = correct.iter().map(|e| Point::new(e[0], e[1])).collect();
        let simplified = visvalingam(&points_ls, &0.0005);
        assert_eq!(simplified, correct_ls);
    }
    #[test]
    fn visvalingam_test_empty_linestring() {
        let vec = Vec::new();
        let compare = Vec::new();
        let simplified = visvalingam(&vec, &1.0);
        assert_eq!(simplified, compare);
    }
    #[test]
    fn visvalingam_test_two_point_linestring() {
        let mut vec = Vec::new();
        vec.push(Point::new(0.0, 0.0));
        vec.push(Point::new(27.8, 0.1));
        let mut compare = Vec::new();
        compare.push(Point::new(0.0, 0.0));
        compare.push(Point::new(27.8, 0.1));
        let simplified = visvalingam(&vec, &1.0);
        assert_eq!(simplified, compare);
    }

    #[test]
    fn multilinestring() {
        // this is the PostGIS example
        let points = vec![
            (5.0, 2.0),
            (3.0, 8.0),
            (6.0, 20.0),
            (7.0, 25.0),
            (10.0, 10.0),
        ];
        let points_ls: Vec<_> = points.iter().map(|e| Point::new(e.0, e.1)).collect();

        let correct = vec![(5.0, 2.0), (7.0, 25.0), (10.0, 10.0)];
        let correct_ls: Vec<_> = correct.iter().map(|e| Point::new(e.0, e.1)).collect();

        let mline = MultiLineString(vec![LineString(points_ls)]);
        assert_eq!(
            mline.simplifyvw(&30.),
            MultiLineString(vec![LineString(correct_ls)])
        );
    }

    #[test]
    fn polygon() {
        let poly = Polygon::new(
            LineString(vec![
                Point::new(0., 0.),
                Point::new(0., 10.),
                Point::new(5., 11.),
                Point::new(10., 10.),
                Point::new(10., 0.),
                Point::new(0., 0.),
            ]),
            vec![],
        );

        let poly2 = poly.simplifyvw(&10.);

        assert_eq!(
            poly2,
            Polygon::new(
                LineString(vec![
                    Point::new(0., 0.),
                    Point::new(0., 10.),
                    Point::new(10., 10.),
                    Point::new(10., 0.),
                    Point::new(0., 0.),
                ]),
                vec![],
            )
        );
    }

    #[test]
    fn multipolygon() {
        let mpoly = MultiPolygon(vec![
            Polygon::new(
                LineString(vec![
                    Point::new(0., 0.),
                    Point::new(0., 10.),
                    Point::new(5., 11.),
                    Point::new(10., 10.),
                    Point::new(10., 0.),
                    Point::new(0., 0.),
                ]),
                vec![]
            ),
        ]);

        let mpoly2 = mpoly.simplifyvw(&10.);

        assert_eq!(
            mpoly2,
            MultiPolygon(vec![
                Polygon::new(
                    LineString(vec![
                        Point::new(0., 0.),
                        Point::new(0., 10.),
                        Point::new(10., 10.),
                        Point::new(10., 0.),
                        Point::new(0., 0.),
                    ]),
                    vec![]
                ),
            ])
        );
    }
}
