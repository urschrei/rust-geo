// Pirzadeh, H. (1999) Computational geometry with the rotating calipers., pp30 – 32
// Available from: http://digitool.library.mcgill.ca/R/-?func=dbin-jump-full&object_id=21623&silo_library=GEN01
// http://web.archive.org/web/20150330010154/http://cgm.cs.mcgill.ca/%7Eorm/rotcal.html
use num_traits::Float;
use types::{Point, Polygon, MultiPolygon, LineString, MultiPoint, MultiLineString};
use std::fmt::Debug;
use std::mem;
// use algorithm::hull_helpers::{swap_remove_to_first, swap_remove_to_last, partition, point_location};
use algorithm::convexhull::ConvexHull;
use algorithm::distance::Distance;
use algorithm::rotate::Rotate;

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
    let mut poly1_hull = poly1.convex_hull().exterior.0.reverse();
    let mut poly2_hull = poly2.convex_hull().exterior.0.reverse();
    println!("Poly 1 Hull: {:?}", poly1.exterior.0);
    println!("Poly 2 Hull {:?}", poly2.exterior.0);
    let (poly1_ymin, poly1_ymax, poly1_xmin, poly1_xmax) =
        (Point::new(6., 1.), Point::new(5., 4.), Point::new(4., 2.), Point::new(7., 3.));
    let (poly2_ymin, poly2_ymax, poly2_xmin, poly2_xmax) =
        (Point::new(11., 1.), Point::new(10., 5.), Point::new(9., 2.), Point::new(12., 3.));

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

    let p_next = Point::new(5., 4.);
    let q_next = Point::new(12., 2.);

    let lower_next = vector_angle(&poly1_ymin,
                                  &Point::new(poly1_xmin.x(), poly1_ymin.y()),
                                  &poly1_ymin,
                                  &p_next);

    let upper_next = vector_angle(&poly2_ymax,
                                  &Point::new(poly2_xmax.x(), poly2_ymax.y()),
                                  &poly2_ymax,
                                  &q_next);

    println!("poly 1 min Y: {:?}", poly1_ymin.y());
    println!("poly 1 min Y, x component: {:?}", poly1_ymin.x());
    println!("poly 2 max Y: {:?}", poly2_ymax.y());
    println!("poly 2 max Y, x component: {:?}", poly2_ymax.x());
    println!("Bottom support (r to l: {:?}", lpoly_1);
    println!("Top support (l to r): {:?}", lpoly_2);
    println!("Minimum distance: {:?}", mindist);
    println!("Lower Angle: {:?}", lower_angle);
    println!("Upper Angle: {:?}", upper_angle);
    println!("Lower Next: {:?}", lower_next);
    println!("Upper Next: {:?}", upper_next);
    println!("Minimum: {:?}", lower_angle.min(upper_angle));

    let rotated = lpoly_1.rotate(lower_angle.min(upper_angle), &poly1_ymin);
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

fn unitvector<T>(slope: T, poly: &Polygon<T>, p: &Point<T>) -> Point<T>
    where T: Float
{
    let clockwise = true;
    let tansq = slope * slope;
    let cossq = T::one() / (T::one() + tansq);
    let sinsq = T::one() - cossq;
    let mut cos = T::zero();
    let mut sin = T::zero();
    // not sure whether this is correct!
    let polysize = poly.exterior.0.len();
    // vertexAt((p.n + T::one())%polysize)
    let pnext_idx = (poly.exterior
                         .0
                         .iter()
                         .position(|&point| point == *p)
                         .unwrap() + 1) % polysize;
    let pnext: Point<T> = poly.exterior.0[pnext_idx];
    // vertexAt((p.n + polysize - T::one())%polysize);
    let pprev_idx = (poly.exterior
                         .0
                         .iter()
                         .position(|&point| point == *p)
                         .unwrap() + polysize - 1) % polysize;
    let pprev: Point<T> = poly.exterior.0[pprev_idx];
    if slope != T::zero() {
        cos = cossq.sqrt();
        sin = sinsq.sqrt();
        if pnext.x() > p.x() {
            if pprev.x() > p.x() {
                if pprev.y() >= p.y() && pnext.y() >= p.y() {
                    if slope > T::zero() {
                        let (mut slprev, mut slnext) = (T::zero(), T::zero());
                        slprev = (pprev.y() - p.y()) / (pprev.x() - p.x());
                        slnext = (pnext.y() - p.y()) / (pnext.x() - p.x());
                        if clockwise && slope <= slprev || !clockwise && slope >= slprev {
                            cos = -cos;
                            sin = -sin;
                        }
                    } else if clockwise {
                        cos = -cos;
                    } else {
                        sin = -sin;
                    }
                } else if pprev.y() <= p.y() && pnext.y() <= p.y() {
                    if slope > T::zero() {
                        if !clockwise {
                            cos = -cos;
                            sin = -sin;
                        } else {
                            let (mut slprev, mut slnext) = (T::zero(), T::zero());
                            slprev = (pprev.y() - p.y()) / (pprev.x() - p.x());
                            slnext = (pnext.y() - p.y()) / (pnext.x() - p.x());
                            if clockwise {
                                if slope <= slprev {
                                    cos = -cos;
                                } else {
                                    sin = -sin;
                                }
                            } else if slope <= slnext {
                                sin = -sin;
                            } else {
                                cos = -cos;
                            }
                        }
                    } else if slope > T::zero() {
                        if !clockwise {
                            cos = -cos;
                            sin = -sin;
                        } else if clockwise {
                            cos = -cos;
                        } else {
                            sin = -sin;
                        }
                    }
                }
            } else if slope < T::zero() {
                //pprev.x() <= p.x()
                sin = -sin;
            }
        } else if pnext.x() < p.x() {
            //pnext.x() <= p.x()
            if pprev.x() < p.x() {
                if (pprev.y() >= p.y()) && (pnext.y() >= p.y()) {
                    if slope > T::zero() {
                        if clockwise {
                            cos = -cos;
                            sin = -sin;
                        } else {
                            let (mut slprev, mut slnext) = (T::zero(), T::zero());
                            slprev = (p.y() - pprev.y()) / (p.x() - pprev.x());
                            slnext = (p.y() - pnext.y()) / (p.x() - pnext.x());
                            if clockwise {
                                if slope <= slprev {
                                    sin = -sin;
                                } else {
                                    cos = -cos;
                                }
                            } else if slope <= slnext {
                                cos = -cos;
                            } else {
                                sin = -sin;
                            }
                        }
                    } else if pprev.y() <= p.y() && pnext.y() <= p.y() {
                        if slope > T::zero() {
                            let (mut slprev, mut slnext) = (T::zero(), T::zero());
                            slprev = (p.y() - pprev.y()) / (p.x() - pprev.x());
                            slnext = (p.y() - pnext.y()) / (p.x() - pnext.x());
                            if slope >= slnext {
                                cos = -cos;
                                sin = -sin;
                            } else if clockwise {
                                sin = -sin;
                            } else {
                                cos = -cos;
                            }
                        } else if slope > T::zero() {
                            if clockwise {
                                cos = -cos;
                                sin = -sin;
                            }
                        } else if clockwise {
                            sin = -sin;
                        } else {
                            cos = -cos;
                        }
                    }
                } else {
                    //pprev.x() >= p.x()
                    cos = -cos;
                    if slope > T::zero() {
                        sin = -sin;
                    }
                }
            } else if pprev.x() > p.x() {
                cos = -cos;
                if slope > T::zero() {
                    sin = -sin;
                }
            } else if slope < T::zero() {
                sin = -sin;
            }
        }
    } else {
        //slope is T::zero()
        sin = T::zero();
        if pnext.x() > p.x() {
            cos = T::one();
        } else if pnext.x() < p.x() {
            cos = -T::one();
        } else if pnext.x() == p.x() {
            if pprev.x() < p.x() {
                cos = T::one();
            } else {
                //pprev > p.x() (can't be equal)
                cos = -T::one();
            }
        }
    }
    Point::new(p.x() + T::from(100).unwrap() * cos,
               p.y() + T::from(100).unwrap() * sin)
}

#[cfg(test)]
mod test {
    use types::Point;
    use super::*;
    #[test]
    fn test_minimum_polygon_distance() {
        let points_raw = vec![(5., 1.), (4., 2.), (4., 3.), (5., 4.), (6., 4.), (7., 3.),
                              (7., 2.), (6., 1.), (5., 1.)];
        let mut points = points_raw
            .iter()
            .map(|e| Point::new(e.0, e.1))
            .collect::<Vec<_>>();
        let poly1 = Polygon::new(LineString(points), vec![]);

        let points_raw_2 = vec![(10., 1.), (9., 2.), (9., 3.), (10., 5.), (11., 4.), (12., 3.),
                                (12., 2.), (11., 1.), (10., 1.)];
        let mut points2 = points_raw_2
            .iter()
            .map(|e| Point::new(e.0, e.1))
            .collect::<Vec<_>>();
        let poly2 = Polygon::new(LineString(points2), vec![]);
        let dist = min_polygon_distance(poly1.convex_hull(), poly2.convex_hull());
        assert_eq!(dist, 3.0);
    }
}
