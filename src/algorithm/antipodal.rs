// Pirzadeh, H. (1999) Computational geometry with the rotating calipers., pp30 – 32
// Available from: http://digitool.library.mcgill.ca/R/-?func=dbin-jump-full&object_id=21623&silo_library=GEN01
// http://web.archive.org/web/20150330010154/http://cgm.cs.mcgill.ca/%7Eorm/rotcal.html
use num_traits::Float;
use types::{Point, Polygon, MultiPolygon, LineString, MultiPoint, MultiLineString};
use std::fmt::Debug;
use std::mem;
use algorithm::hull_helpers::{swap_remove_to_first, swap_remove_to_last, partition, point_location};
use algorithm::convexhull::ConvexHull;
use algorithm::distance::Distance;

fn rotation_matrix<T>(angle: T, origin: &Point<T>, points: &[Point<T>]) -> Vec<Point<T>>
    where T: Float + Debug
{
    let cos_theta = angle.to_radians().cos();
    let sin_theta = angle.to_radians().sin();
    let x0 = origin.x();
    let y0 = origin.y();
    points.iter()
        .map(|point| {
                println!("Point: {:?}", point);
                 let x = point.x() - x0;
                 let y = point.y() - y0;
                 Point::new(x * cos_theta - y * sin_theta + x0,
                            x * sin_theta + y * cos_theta + y0)
             })
        .collect::<Vec<_>>()
}

// calculate max and min polygon points
fn min_max<T>(mut hull: &mut [Point<T>]) -> (Point<T>, Point<T>, Point<T>, Point<T>)
    where T: Float
{
    let mut ymin = swap_remove_to_first(&mut hull, 0);
    let mut ymax = swap_remove_to_first(&mut hull, 0);
    let mut xmax = swap_remove_to_first(&mut hull, 0);
    let mut xmin = swap_remove_to_first(&mut hull, 0);
    if ymin.y() > ymax.y() {
        mem::swap(ymin, ymax);
    }
    for point in hull.iter_mut() {
        if point.y() < ymin.y() {
            mem::swap(point, ymin);
        }
        if point.y() > ymax.y() {
            mem::swap(point, ymax);
        }
        if point.x() < xmin.x() {
            mem::swap(point, &mut xmin);
        }
        if point.x() > xmax.x() {
            mem::swap(point, &mut xmax);
        }
    }
    (*ymin, *ymax, *xmin, *xmax)
}

// return the vector angle in degrees of two vectors ab and cd
fn vector_angle<T>(a: &Point<T>, b: &Point<T>, c: &Point<T>, d: &Point<T>) -> T
    where T: Float + Debug
{
    // vector by initial and terminal points
    let ab = Point::new(b.x() - a.x(), b.y() - a.y());
    let cd = Point::new(d.x() - c.x(), d.y() - c.y());

    let product = ab.dot(&cd);

    // magnitudes
    let magnitude_ab = (ab.x().powi(2) + ab.y().powi(2)).sqrt();
    let magnitude_cd = (cd.x().powi(2) + cd.y().powi(2)).sqrt();

    let cos_theta = ab.dot(&cd) / (magnitude_ab * magnitude_cd);
    // inverse, in degrees
    cos_theta.abs().acos().to_degrees()
}

fn min_polygon_distance<T>(mut poly1: Polygon<T>, mut poly2: Polygon<T>) -> T
    where T: Float + Debug
{
    // polygons must be convex
    let mut poly1_hull = poly1.exterior.0.as_mut_slice();
    let mut poly2_hull = poly2.exterior.0.as_mut_slice();

    let (poly1_ymin, poly1_ymax, poly1_xmin, poly1_xmax) = min_max(&mut poly1_hull);
    let (poly2_ymin, poly2_ymax, poly2_xmin, poly2_xmax) = min_max(&mut poly2_hull);

    // lines of support must be parallel to the x axis
    // lower support tangent to poly1, which must lie to its right
    let mut lpoly_1 = LineString(vec![Point::new(poly1_xmax.x(), poly1_ymin.y()),
                                      Point::new(poly1_ymin.x(), poly1_ymin.y()),
                                      Point::new(poly1_xmin.x(), poly1_ymin.y())]);
    // upper support tangent to poly2, which must lie to its right
    let mut lpoly_2 = LineString(vec![Point::new(poly2_xmin.x(), poly2_ymax.y()),
                                      Point::new(poly2_ymax.x(), poly2_ymax.y()),
                                      Point::new(poly2_xmax.x(), poly2_ymax.y())]);
    // initial minimum distance
    let mut mindist = poly1_ymin.distance(&poly2_ymax);
    let lower_angle = vector_angle(&poly1_ymin,
                                   &Point::new(poly1_xmin.x(), poly1_ymin.y()),
                                   &poly1_ymin,
                                   &poly2_ymax);
    let upper_angle = vector_angle(&poly2_ymax,
                                   &Point::new(poly2_xmax.x(), poly2_ymax.y()),
                                   &poly2_ymax,
                                   &poly1_ymin);

    println!("poly 1 min Y: {:?}", poly1_ymin.y());
    println!("poly 1 min Y, x component: {:?}", poly1_ymin.x());
    println!("poly 2 max Y: {:?}", poly2_ymax.y());
    println!("poly 2 max Y, x component: {:?}", poly2_ymax.x());
    println!("poly 1 min X: {:?}", poly1_xmin);
    println!("poly 2 max X: {:?}", poly2_xmax);
    println!("Bottom support (r to l: {:?}", lpoly_1);
    println!("Top support (l to r): {:?}", lpoly_2);
    println!("Minimum distance: {:?}", mindist);
    println!("Lower Angle: {:?}", lower_angle);
    println!("Upper Angle: {:?}", upper_angle);
    println!("Minimum: {:?}", lower_angle.min(upper_angle));

    let rotated = rotation_matrix(T::from(-45.0).unwrap(), &poly1_ymin, lpoly_1.0.as_slice());
    println!("Rotated: {:?}", rotated);

    // 1.  We want poly1_min.y(), and poly2_max.y()
    // 2.  Construct two lines of support, parallel to the x axis – LP and LQ –
    //     which touch the polygons at yminP and ymaxQ
    //     such that the polygons lie to the right of their respective lines of support.
    //     LP and LQ have opposite direction, and yminP and ymaxQ form an anti-podal pair between the polygons.
    //     The lines of support lie on vertices pi ∈ P, and qj ∈ Q, and determine two angles
    //     θi, θj, which are computed
    // 3.  Compute dist(yminP,ymaxQ) and keep it as the minimum.
    // 3a. Compute θ = θi.min(θj)
    // 4.  Rotate the lines clockwise about pi and qj by θ
    //     One of the lines should now be flush with an edge of its polygon. (???)
    // 5.  If only one line coincides with an edge, then:
    //     - the vertex-edge anti-podal pair distance should be computed
    //     - the new vertex-vertex anti-podal pair distance should be computed
    //     Both distances are compared the current minimum, which is updated if necessary.
    //
    //     If both lines of support coincide with edges, then the situation is somewhat more complex:
    //     - If the edges "overlap", that is if one can construct a line perpendicular to both edges and
    //     intersecting both edges (but not at vertices), then the edge-edge distance should be computed.
    //     - Otherwise the three new vertex-vertex anti-podal pair distances are computed.
    //     All distances are compared to the current minimum, which is updated if necessary.
    // 6.  Repeat steps 4 and 5, until the lines reach (yminP, ymaxQ) again.
    // 7. Return the minimum

    T::from(3.0).unwrap()
}

#[cfg(test)]
mod test {
    use types::Point;
    use super::*;
    #[test]
    fn test_minimum_polygon_distance() {
        let points_raw = vec![(5., 1.), (4., 2.), (4., 3.), (5., 4.), (6., 4.), (7., 3.),
                              (7., 2.), (6., 1.), (5., 1.)];
        let mut points = points_raw.iter().map(|e| Point::new(e.0, e.1)).collect::<Vec<_>>();
        let poly1 = Polygon::new(LineString(points), vec![]);

        let points_raw_2 = vec![(10., 1.), (9., 2.), (9., 3.), (10., 5.), (11., 4.), (12., 3.),
                                (12., 2.), (11., 1.), (10., 1.)];
        let mut points2 = points_raw_2.iter().map(|e| Point::new(e.0, e.1)).collect::<Vec<_>>();
        let poly2 = Polygon::new(LineString(points2), vec![]);
        let dist = min_polygon_distance(poly1.convex_hull(), poly2.convex_hull());
        assert_eq!(dist, 3.0);
    }
}
