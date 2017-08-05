use num_traits::Float;
use types::{Point, LineString, Polygon};
use algorithm::distance::Distance;
use algorithm::simplify::Simplify;

fn radial_distance<T>(points: &[Point<T>], epsilon: &T) -> Vec<Point<T>>
where
    T: Float
{
    if points.is_empty() {
        return points.to_vec();
    }
    let last = points.last().unwrap();
    let mut prev_point = points[0];
    let mut new_points = vec![prev_point];
    for &point in points.iter().skip(1) {
        if point.distance(&prev_point).powi(2) > epsilon.powi(2) {
            new_points.push(point);
            prev_point = point;
        }
    }
    if prev_point != *last {
        new_points.push(*last);
    }
    new_points
}

// Radial-distance algorithm
pub trait SimplifyRadial<T, Epsilon = T> {
    /// Returns the simplified representation of a Geometry, using a radial-distance pre-processing step  
    /// This reduces the overall quality of representation, but is faster.
    ///
    /// ```
    /// use geo::{Point, LineString};
    /// use geo::algorithm::simplify_radial::{SimplifyRadial};
    ///
    /// let mut vec = Vec::new();
    /// vec.push(Point::new(0.0, 0.0));
    /// vec.push(Point::new(5.0, 4.0));
    /// vec.push(Point::new(11.0, 5.5));
    /// vec.push(Point::new(17.3, 3.2));
    /// vec.push(Point::new(27.8, 0.1));
    /// let linestring = LineString(vec);
    /// let mut compare = Vec::new();
    /// compare.push(Point::new(0.0, 0.0));
    /// compare.push(Point::new(5.0, 4.0));
    /// compare.push(Point::new(11.0, 5.5));
    /// compare.push(Point::new(27.8, 0.1));
    /// let ls_compare = LineString(compare);
    /// let simplified = linestring.simplify_radial(&1.0);
    /// assert_eq!(simplified, ls_compare)
    /// ```
    fn simplify_radial(&self, epsilon: &T) -> Self
    where
        T: Float;
}

impl<T> SimplifyRadial<T> for LineString<T>
where
    T: Float
{
    fn simplify_radial(&self, epsilon: &T) -> LineString<T> {
        LineString(radial_distance(&self.0, epsilon)).simplify(&epsilon)
    }
}

impl<T> SimplifyRadial<T> for Polygon<T>
where
    T: Float,
{
    fn simplify_radial(&self, epsilon: &T) -> Polygon<T> {
        Polygon::new(
            self.exterior.simplify_radial(epsilon).simplify(epsilon),
            self.interiors
                .iter()
                .map(|ls| ls.simplify_radial(epsilon).simplify(epsilon))
                .collect(),
        )
    }
}

#[cfg(test)]
mod test {
    use types::{Point, LineString, Polygon};
    use algorithm::simplify::Simplify;
    use algorithm::simplify_radial::SimplifyRadial;
    use super::radial_distance;

    #[test]
    fn test_radial() {
        // test values stolen from simplify.js
        let points = vec![
            Point::new(224.55,250.15), Point::new(226.91,244.19), Point::new(233.31,241.45), Point::new(234.98,236.06),
            Point::new(244.21,232.76), Point::new(262.59,215.31), Point::new(267.76,213.81), Point::new(273.57,201.84),
            Point::new(273.12,192.16), Point::new(277.62,189.03), Point::new(280.36,181.41), Point::new(286.51,177.74),
            Point::new(292.41,159.37), Point::new(296.91,155.64), Point::new(314.95,151.37), Point::new(319.75,145.16),
            Point::new(330.33,137.57), Point::new(341.48,139.96), Point::new(369.98,137.89), Point::new(387.39,142.51),
            Point::new(391.28,139.39), Point::new(409.52,141.14), Point::new(414.82,139.75), Point::new(427.72,127.30),
            Point::new(439.60,119.74), Point::new(474.93,107.87), Point::new(486.51,106.75), Point::new(489.20,109.45),
            Point::new(493.79,108.63), Point::new(504.74,119.66), Point::new(512.96,122.35), Point::new(518.63,120.89),
            Point::new(524.09,126.88), Point::new(529.57,127.86), Point::new(534.21,140.93), Point::new(539.27,147.24),
            Point::new(567.69,148.91), Point::new(575.25,157.26), Point::new(580.62,158.15), Point::new(601.53,156.85),
            Point::new(617.74,159.86), Point::new(622.00,167.04), Point::new(629.55,194.60), Point::new(638.90,195.61),
            Point::new(641.26,200.81), Point::new(651.77,204.56), Point::new(671.55,222.55), Point::new(683.68,217.45),
            Point::new(695.25,219.15), Point::new(700.64,217.98), Point::new(703.12,214.36), Point::new(712.26,215.87),
            Point::new(721.49,212.81), Point::new(727.81,213.36), Point::new(729.98,208.73), Point::new(735.32,208.20),
            Point::new(739.94,204.77), Point::new(769.98,208.42), Point::new(779.60,216.87), Point::new(784.20,218.16),
            Point::new(800.24,214.62), Point::new(810.53,219.73), Point::new(817.19,226.82), Point::new(820.77,236.17),
            Point::new(827.23,236.16), Point::new(829.89,239.89), Point::new(851.00,248.94), Point::new(859.88,255.49),
            Point::new(865.21,268.53), Point::new(857.95,280.30), Point::new(865.48,291.45), Point::new(866.81,298.66),
            Point::new(864.68,302.71), Point::new(867.79,306.17), Point::new(859.87,311.37), Point::new(860.08,314.35),
            Point::new(858.29,314.94), Point::new(858.10,327.60), Point::new(854.54,335.40), Point::new(860.92,343.00),
            Point::new(856.43,350.15), Point::new(851.42,352.96), Point::new(849.84,359.59), Point::new(854.56,365.53),
            Point::new(849.74,370.38), Point::new(844.09,371.89), Point::new(844.75,380.44), Point::new(841.52,383.67),
            Point::new(839.57,390.40), Point::new(845.59,399.05), Point::new(848.40,407.55), Point::new(843.71,411.30),
            Point::new(844.09,419.88), Point::new(839.51,432.76), Point::new(841.33,441.04), Point::new(847.62,449.22),
            Point::new(847.16,458.44), Point::new(851.38,462.79), Point::new(853.97,471.15), Point::new(866.36,480.77)
        ];
        let correct = vec![
            Point::new(224.55,250.15), Point::new(267.76,213.81), Point::new(296.91,155.64), Point::new(330.33,137.57),
            Point::new(409.52,141.14), Point::new(439.60,119.74), Point::new(486.51,106.75), Point::new(529.57,127.86),
            Point::new(539.27,147.24), Point::new(617.74,159.86), Point::new(629.55,194.60), Point::new(671.55,222.55),
            Point::new(727.81,213.36), Point::new(739.94,204.77), Point::new(769.98,208.42), Point::new(779.60,216.87),
            Point::new(800.24,214.62), Point::new(820.77,236.17), Point::new(859.88,255.49), Point::new(865.21,268.53),
            Point::new(857.95,280.30), Point::new(867.79,306.17), Point::new(859.87,311.37), Point::new(854.54,335.40),
            Point::new(860.92,343.00), Point::new(849.84,359.59), Point::new(854.56,365.53), Point::new(844.09,371.89),
            Point::new(839.57,390.40), Point::new(848.40,407.55), Point::new(839.51,432.76), Point::new(853.97,471.15),
            Point::new(866.36,480.77)
        ];
        let simplified = LineString(radial_distance(&points, &5.0)).simplify(&5.0);
        assert_eq!(simplified.0, correct);
    }
    #[test]
    fn radial_distance_test_empty_linestring() {
        let vec = Vec::new();
        let compare = Vec::new();
        let simplified = radial_distance(&vec, &5.0);
        assert_eq!(simplified, compare);
    }
    #[test]
    fn radial_distance_test_two_point_linestring() {
        let mut vec = Vec::new();
        vec.push(Point::new(0.0, 0.0));
        vec.push(Point::new(27.8, 0.1));
        let mut compare = Vec::new();
        compare.push(Point::new(0.0, 0.0));
        compare.push(Point::new(27.8, 0.1));
        let simplified = radial_distance(&vec, &5.0);
        assert_eq!(simplified, compare);
    }
    #[test]
    fn polygon_simplification_test() {
        let ls1 = LineString(vec![
            Point::new(5.0, 1.0),
            Point::new(4.0, 2.0),
            Point::new(4.0, 3.0),
            Point::new(5.0, 4.0),
            Point::new(6.0, 4.0),
            Point::new(7.0, 3.0),
            Point::new(7.0, 2.0),
            Point::new(6.0, 1.0),
            Point::new(5.0, 1.0),
        ]);

        let ls2 = LineString(vec![
            Point::new(5.0, 1.3),
            Point::new(5.5, 2.0),
            Point::new(6.0, 1.3),
            Point::new(5.0, 1.3),
        ]);

        let correct_outside = vec![
            (5.0, 1.0),
            (4.0, 3.0),
            (6.0, 4.0),
            (7.0, 2.0),
            (6.0, 1.0),
            (5.0, 1.0),
        ].iter()
            .map(|e| Point::new(e.0, e.1))
            .collect::<Vec<_>>();

        let poly1 = Polygon::new(ls1, vec![ls2]);
        let simplified = poly1.simplify_radial(&0.45);

        assert_eq!(simplified.exterior.0, correct_outside);
    }
}
