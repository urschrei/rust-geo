// Pirzadeh, H. (1999) Computational geometry with the rotating calipers., pp30 – 32
// Available from: http://digitool.library.mcgill.ca/R/-?func=dbin-jump-full&object_id=21623&silo_library=GEN01
use num_traits::Float;
use types::{Point, Polygon, MultiPolygon, LineString, MultiPoint, MultiLineString};
use std::fmt::Debug;
use std::mem;
use algorithm::hull_helpers::{
    swap_remove_to_first,
    partition,
    point_location
};

fn antipodal<T>(mut points: &mut [Point<T>]) -> Vec<(Point<T>, Point<T>)>
    where T: Float + Debug
{
    // can't build a hull from fewer than four points
    if points.len() < 4 {
        // return 0 as usize
        return vec![(Point::new(T::zero(), T::zero()), Point::new(T::zero(), T::zero()))];
    }
    let mut min = swap_remove_to_first(&mut points, 0);
    let mut max = swap_remove_to_first(&mut points, 0);
    if min.x() > max.x() {
        mem::swap(min, max);
    }
    for point in points.iter_mut() {
        if point.x() < min.x() {
            mem::swap(point, min);
        }
        if point.x() > max.x() {
            mem::swap(point, max);
        }
    }
    let last_upper = partition(&mut points, |p| point_location(max, min, p));
    let upper = points[..last_upper].to_vec();
    let last_lower = partition(&mut points, |p| point_location(min, max, p));
    let lower = points[..last_lower].to_vec();
    let mut antipodal: Vec<(Point<T>, Point<T>)> = vec![];
    let mut i = 0;
    let mut j = lower.len() - 1;
    while i < upper.len() - 1 || j > 0 {
        antipodal.push((upper[i], lower[j]));
        if i == upper.len() - 1 {
            // we've walked all the way along the upper hull
            j -= 1;
        } else if j == 0 {
            // we've walked all the way along the lower hull
            i += 1;
            // There are points remaining, so compare the slopes of next hull edges
            // Could we have divide-by-0 errors here?
        } else if (upper[i + 1].y() - upper[i].y()) * (lower[j].x() - lower[j - 1].x()) >
                  (upper[i + 1].x() - upper[i].x()) * (lower[j].y() - lower[j - 1].y()) {
            i += 1;
        } else {
            j -= 1;
        }
    }
    antipodal
}

#[cfg(test)]
mod test {
    use types::Point;
    use super::*;
    #[test]
    fn test_antipodes() {
        let points_raw = vec![(5., 1.), (4., 2.), (4., 3.), (5., 4.), (6., 4.), (7., 3.),
                              (7., 2.), (6., 1.), (5., 1.)];
        let mut points = points_raw.iter().map(|e| Point::new(e.0, e.1)).collect::<Vec<_>>();
        // currently have no idea whether this result is correct
        let antipodes = antipodal(points.as_mut_slice());
        let correct = vec![(Point::new(5.0, 1.0), Point::new(6.0, 4.0)),
                           (Point::new(5.0, 1.0), Point::new(5.0, 4.0)),
                           (Point::new(5.0, 1.0), Point::new(4.0, 3.0)),
                           (Point::new(6.0, 1.0), Point::new(4.0, 3.0)),
                           (Point::new(5.0, 1.0), Point::new(4.0, 3.0))];
        assert_eq!(antipodes, correct);
    }
}
