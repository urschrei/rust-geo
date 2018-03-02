use num_traits::{Float, FromPrimitive};
use types::Point;

/// Returns a new Point using the distance to the existing Point and a bearing for the direction

pub trait HaversineDestination<T: Float> {
    /// Returns a new Point using distance to the existing Point and a bearing for the direction
    ///
    /// ```
    /// use geo::Point;
    /// use geo::algorithm::haversine_destination::HaversineDestination;
    ///
    /// let p_1 = Point::<f64>::new(9.177789688110352, 48.776781529534965);
    /// let p_2 = p_1.haversine_destination(45., 10000.);
    /// assert_eq!(p_2, Point::<f64>::new(9.274410083250379, 48.84033282787534))
    /// ```
    fn haversine_destination(&self, bearing: T, distance: T) -> Point<T>;
}

impl<T> HaversineDestination<T> for Point<T>
where
    T: Float + FromPrimitive,
{
    fn haversine_destination(&self, bearing: T, distance: T) -> Point<T> {
        let center_lng = self.x().to_radians();
        let center_lat = self.y().to_radians();
        let bearing_rad = bearing.to_radians();

        // WGS84 equatorial radius is 6378137.0
        let rad = distance / T::from(6371000.0).unwrap();

        let lat = {
            center_lat.sin() * rad.cos() + center_lat.cos() * rad.sin() * bearing_rad.cos()
        }.asin();
        let lng = { bearing_rad.sin() * rad.sin() * center_lat.cos() }
            .atan2(rad.cos() - center_lat.sin() * lat.sin()) + center_lng;

        Point::new(lng.to_degrees(), lat.to_degrees())
    }
}

pub trait HaversineBearing<T: Float> {
    /// Returns a new bearing, given origin and destination points
    ///
    /// ```
    /// use geo::Point;
    /// use geo::algorithm::haversine_destination::HaversineBearing;
    ///
    /// let orig = Point::<f64>::new(0.00000002, 51.23565754);
    /// let dest = Point::<f64>::new(0.00000012, 51.23565764);
    /// let bearing = haversine_bearing(&orig, &dest);
    /// assert_eq!(bearing, 10.0)
    /// ```
    fn haversine_bearing(&self, origin: &Point<T>, destination: &Point<T>) -> T;
}

impl<T> HaversineBearing<T> for Point<T>
where
    T: Float + FromPrimitive,
{
    fn haversine_bearing(&self, origin: &Point<T>, destination: &Point<T>) -> T {
        // https://gist.github.com/geografa/1366401
        let (lon1, lat1) = (origin.x().to_radians(), origin.y().to_radians());
        let (lon2, lat2) = (destination.x().to_radians(), destination.y().to_radians());
        let dlon = lon2 - lon1;
        let bearing = T::atan2(
            T::sin(dlon) * T::cos(lat2),
            T::cos(lat1) * T::sin(lat2) - T::sin(lat1) * T::cos(lat2) * T::cos(dlon),
        );
        bearing.to_degrees()
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use algorithm::haversine_distance::{HaversineDistance};
    use num_traits::pow;

    #[test]
    fn returns_a_new_point() {
        let p_1 = Point::<f64>::new(9.177789688110352, 48.776781529534965);
        let p_2 = p_1.haversine_destination(45., 10000.);
        assert_eq!(p_2, Point::<f64>::new(9.274410083250379, 48.84033282787534));
        let distance = p_1.haversine_distance(&p_2);
        assert_relative_eq!(distance, 10000., epsilon = 1.0e-6)
    }

    #[test]
    fn direct_and_indirect_destinations_are_close() {
        let p_1 = Point::<f64>::new(9.177789688110352, 48.776781529534965);
        let p_2 = p_1.haversine_destination(45., 10000.);
        let square_edge = { pow(10000., 2) / 2. }.sqrt();
        let p_3 = p_1.haversine_destination(0., square_edge);
        let p_4 = p_3.haversine_destination(90., square_edge);
        assert_relative_eq!(p_4.x(), p_2.x(), epsilon = 1.0e-6);
        assert_relative_eq!(p_4.y(), p_2.y(), epsilon = 1.0e-6);
    }
}
