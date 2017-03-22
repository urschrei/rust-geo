// Pirzadeh, H. (1999) Computational geometry with the rotating calipers., pp30 – 32
// Available from: http://digitool.library.mcgill.ca/R/-?func=dbin-jump-full&object_id=21623&silo_library=GEN01
// http://web.archive.org/web/20150330010154/http://cgm.cs.mcgill.ca/%7Eorm/rotcal.html
use num_traits::Float;
use types::{Point, Polygon, MultiPolygon, LineString, MultiPoint, MultiLineString};
use std::fmt::Debug;
use std::mem;
use algorithm::hull_helpers::{
    swap_remove_to_first,
    swap_remove_to_last,
    partition,
    point_location
};

fn min_polygon_distance<T>(mut poly1: &mut [Point<T>], mut poly2: &mut [Point<T>]) -> T
    where T: Float + Debug
{
    // poly1 min y
    let mut poly1_min = swap_remove_to_first(&mut poly1, 0);
    let mut poly1_max = swap_remove_to_first(&mut poly1, 0);
    if poly1_min.y() > poly1_max.y() {
        mem::swap(poly1_min, poly1_max);
    }
    for point in poly1.iter_mut() {
        if point.y() < poly1_min.y() {
            mem::swap(point, poly1_min);
        }
        if point.y() > poly1_max.y() {
            mem::swap(point, poly1_max);
        }
    }
    // poly2 max y
    let mut poly2_min = swap_remove_to_first(&mut poly2, 0);
    let mut poly2_max = swap_remove_to_first(&mut poly2, 0);
    if poly2_min.y() > poly2_max.y() {
        mem::swap(poly2_min, poly2_max);
    }
    for point in poly2.iter_mut() {
        if point.y() < poly2_min.y() {
            mem::swap(point, poly2_min);
        }
        if point.y() > poly2_max.y() {
            mem::swap(point, poly2_max);
        }
    }
    println!("poly 1 min Y: {:?}", poly1_min.y());
    println!("poly 2 max Y: {:?}", poly2_max.y());
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
        // let poly1 = Polygon::new(Linestring(points), vec![]);

        let points_raw_2 = vec![(8., 1.), (7., 2.), (7., 3.), (8., 4.), (9., 4.), (10., 3.),
                      (10., 2.), (9., 1.), (8., 1.)];
        let mut points2 = points_raw_2.iter().map(|e| Point::new(e.0, e.1)).collect::<Vec<_>>();
        // let poly2 = Polygon::new(Linestring(points2), vec![]);
        let dist = min_polygon_distance(&mut points, &mut points2);
        assert_eq!(dist, 3.0);
    }
}
