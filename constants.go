package geodesy

import "math"

const (
	// WGS-84 ellipsoid constants
	// https://skybrary.aero/sites/default/files/bookshelf/5854.pdf
	wgs84A  = 6378137.0 // semi-major axis (m)
	wgs84F  = 1.0 / 298.257223563
	wgs84B  = wgs84A * (1 - wgs84F)
	wgs84E2 = (wgs84A*wgs84A - wgs84B*wgs84B) / (wgs84A * wgs84A)

	// EarthRadius is the mean Earth radius in meters.
	EarthRadius = 6371000.0
)

// deg2rad converts degrees to radians.
func deg2rad(deg float64) float64 {
	return deg * math.Pi / 180.0
}
