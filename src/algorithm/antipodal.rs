use num_traits::Float;
use types::{Point, Polygon, MultiPolygon, LineString, MultiPoint, MultiLineString};
use std::mem;
use std::fmt::Debug;

fn swap_remove_to_first<'a, T>(slice: &mut &'a mut [T], idx: usize) -> &'a mut T {
    let tmp = mem::replace(slice, &mut []);
    tmp.swap(0, idx);
    let (h, t) = tmp.split_first_mut().unwrap();
    *slice = t;
    h
}

fn swap_remove_to_last<'a, T>(slice: &mut &'a mut [T], idx: usize) -> &'a mut T {
    let tmp = mem::replace(slice, &mut []);
    let len = tmp.len();
    tmp.swap(len - 1, idx);
    let (h, t) = tmp.split_last_mut().unwrap();
    *slice = t;
    h
}

// slice[..result] have pred(e) == true, slice[result..] have pred(e) == false
fn partition<T, F: FnMut(&T) -> bool>(mut slice: &mut [T], mut pred: F) -> usize {
    let mut i = 0;
    loop {
        let test = match slice.first() {
            Some(e) => pred(e),
            None => break,
        };
        if test {
            swap_remove_to_first(&mut slice, 0);
            i += 1;
        } else {
            swap_remove_to_last(&mut slice, 0);
        }
    }
    i
}

// Determine whether a point lies on one side of a line segment, or the other.
// The cross product v x w of two vectors v and w is a vector whose length is
// |v||w|sin φ, (where |v| is the length of v and φ is the angle between the vectors),
// and which is orthogonal (perpendicular) to both v and w.  Since there are two
// such possible vectors, the definition arbitrarily selects the one that matches
// the direction in which a screw would move if rotated from v to w

// Mathematically, if the coordinates of vectors v and w are (vx, vy) and (wx, wy)
// respectively, the cross product will be (vxwy - vywx). If a segment is
// defined by points A B and, we wish to check on which side of AB a third point C falls,
// we can compute the cross product AB x AC and check its sign:
// If it's negative, it will be on the "right" side of AB
// (when standing on A and looking towards B). If positive, it will be on the left side
fn cross_prod<T>(p_a: &Point<T>, p_b: &Point<T>, p_c: &Point<T>) -> T
    where T: Float
{
    (p_b.x() - p_a.x()) * (p_c.y() - p_a.y()) - (p_b.y() - p_a.y()) * (p_c.x() - p_a.x())
}
fn point_location<T>(p_a: &Point<T>, p_b: &Point<T>, p_c: &Point<T>) -> bool
    where T: Float
{
    cross_prod(p_a, p_b, p_c) > T::zero()
}

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
            j = j - 1;
        } else if j == 0 {
            // we've walked all the way along the lower hull
            i = i + 1;
            // There are points remaining, so compare the slopes of next hull edges
            // Could we have divide-by-0 errors here?
        } else if (upper[i + 1].y() - upper[i].y()) * (lower[j].x() - lower[j - 1].x()) >
                  (upper[i + 1].x() - upper[i].x()) * (lower[j].y() - lower[j - 1].y()) {
            i = i + 1;
        } else {
            j = j - 1;
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
