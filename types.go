package geodesy

// LLA represents geodetic coordinates on the WGS-84 ellipsoid.
//
// Fields:
//   - Lat: latitude in degrees, positive north (-90 to +90)
//   - Lon: longitude in degrees, positive east (-180 to +180)
//   - Alt: altitude in meters above the WGS-84 ellipsoid (not MSL)
type LLA struct {
	Lat, Lon, Alt float64
}

// Vec3 represents a 3D Cartesian vector with X, Y, Z components.
type Vec3 struct {
	X, Y, Z float64
}

// ECEF represents Earth-Centered Earth-Fixed Cartesian coordinates.
// This is a global coordinate system with origin at Earth's center of mass.
//
// Axes:
//   - X: passes through the equator at the prime meridian (0N, 0E)
//   - Y: passes through the equator at 90E longitude (0N, 90E)
//   - Z: passes through the North Pole
//
// All values are in meters.
type ECEF Vec3

// ENU represents local East-North-Up tangent plane coordinates.
// This is a local coordinate system relative to a reference point on Earth's surface.
//
// Axes:
//   - X (East): positive towards geographic east
//   - Y (North): positive towards geographic north
//   - Z (Up): positive away from Earth's center (perpendicular to ellipsoid)
//
// All values are in meters, representing displacement from the reference origin.
type ENU Vec3
