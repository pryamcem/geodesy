package geodesy

import "math"

// LLAtoECEF converts geodetic coordinates to Earth-Centered Earth-Fixed (ECEF)
// Cartesian coordinates.
// https://en.wikipedia.org/wiki/Geographic_coordinate_conversion#From_geodetic_to_ECEF_coordinates
func LLAtoECEF(p LLA) ECEF {
	lat := deg2rad(p.Lat)
	lon := deg2rad(p.Lon)

	sinLat := math.Sin(lat)
	cosLat := math.Cos(lat)
	sinLon := math.Sin(lon)
	cosLon := math.Cos(lon)

	// Radius of curvature in the prime vertical
	N := wgs84A / math.Sqrt(1.0-wgs84E2*sinLat*sinLat)

	return ECEF{
		X: (N + p.Alt) * cosLat * cosLon,
		Y: (N + p.Alt) * cosLat * sinLon,
		Z: (N*(1.0-wgs84E2) + p.Alt) * sinLat,
	}
}

// ECEFtoENU converts ECEF coordinates to local East-North-Up (ENU) coordinates
// relative to a reference point.
// https://en.wikipedia.org/wiki/Geographic_coordinate_conversion#From_ECEF_to_ENU
func ECEFtoENU(p, origin ECEF, originLLA LLA) ENU {
	lat := deg2rad(originLLA.Lat)
	lon := deg2rad(originLLA.Lon)

	dx := p.X - origin.X
	dy := p.Y - origin.Y
	dz := p.Z - origin.Z

	sinLat := math.Sin(lat)
	cosLat := math.Cos(lat)
	sinLon := math.Sin(lon)
	cosLon := math.Cos(lon)

	return ENU{
		X: -sinLon*dx + cosLon*dy,
		Y: -sinLat*cosLon*dx - sinLat*sinLon*dy + cosLat*dz,
		Z: cosLat*cosLon*dx + cosLat*sinLon*dy + sinLat*dz,
	}
}

// LLAtoENU converts geodetic coordinates directly to local ENU coordinates
// relative to a reference point.
func LLAtoENU(p, origin LLA) ENU {
	ecef := LLAtoECEF(p)
	originECEF := LLAtoECEF(origin)

	return ECEFtoENU(ecef, originECEF, origin)
}

// ENUtoECEF converts local ENU coordinates back to ECEF coordinates
// relative to a reference point.
func ENUtoECEF(p ENU, origin ECEF, originLLA LLA) ECEF {
	lat := deg2rad(originLLA.Lat)
	lon := deg2rad(originLLA.Lon)

	sinLat := math.Sin(lat)
	cosLat := math.Cos(lat)
	sinLon := math.Sin(lon)
	cosLon := math.Cos(lon)

	// Inverse rotation matrix (transpose of ECEFtoENU rotation)
	dx := -sinLon*p.X - sinLat*cosLon*p.Y + cosLat*cosLon*p.Z
	dy := cosLon*p.X - sinLat*sinLon*p.Y + cosLat*sinLon*p.Z
	dz := cosLat*p.Y + sinLat*p.Z

	return ECEF{
		X: origin.X + dx,
		Y: origin.Y + dy,
		Z: origin.Z + dz,
	}
}

// ECEFtoLLA converts ECEF coordinates to geodetic coordinates.
// Uses iterative method for accurate results.
func ECEFtoLLA(p ECEF) LLA {
	// Distance from Z-axis
	rxy := math.Sqrt(p.X*p.X + p.Y*p.Y)

	// Longitude
	lon := math.Atan2(p.Y, p.X)

	// Iterative solution for latitude
	lat := math.Atan2(p.Z, rxy)
	for i := 0; i < 10; i++ {
		sinLat := math.Sin(lat)
		N := wgs84A / math.Sqrt(1.0-wgs84E2*sinLat*sinLat)
		lat = math.Atan2(p.Z+wgs84E2*N*sinLat, rxy)
	}

	// Altitude
	sinLat := math.Sin(lat)
	N := wgs84A / math.Sqrt(1.0-wgs84E2*sinLat*sinLat)
	alt := rxy/math.Cos(lat) - N

	return LLA{
		Lat: lat * 180 / math.Pi,
		Lon: lon * 180 / math.Pi,
		Alt: alt,
	}
}

// ENUtoLLA converts local ENU coordinates back to geodetic coordinates
// relative to a reference point.
func ENUtoLLA(p ENU, origin LLA) LLA {
	originECEF := LLAtoECEF(origin)
	ecef := ENUtoECEF(p, originECEF, origin)
	return ECEFtoLLA(ecef)
}

// ENUtoNED converts ENU coordinates to NED coordinates.
// The frames share the same origin; the transformation is a permutation with sign flip on Z.
func ENUtoNED(p ENU) NED {
	return NED{X: p.Y, Y: p.X, Z: -p.Z}
}

// NEDtoENU converts NED coordinates to ENU coordinates.
func NEDtoENU(p NED) ENU {
	return ENU{X: p.Y, Y: p.X, Z: -p.Z}
}

// LLAtoNED converts geodetic coordinates directly to local NED coordinates
// relative to a reference point.
func LLAtoNED(p, origin LLA) NED {
	return ENUtoNED(LLAtoENU(p, origin))
}

// NEDtoLLA converts local NED coordinates back to geodetic coordinates
// relative to a reference point.
func NEDtoLLA(p NED, origin LLA) LLA {
	return ENUtoLLA(NEDtoENU(p), origin)
}

// ECEFtoNED converts ECEF coordinates to local NED coordinates
// relative to a reference point.
func ECEFtoNED(p, origin ECEF, originLLA LLA) NED {
	return ENUtoNED(ECEFtoENU(p, origin, originLLA))
}

// NEDtoECEF converts local NED coordinates back to ECEF coordinates
// relative to a reference point.
func NEDtoECEF(p NED, origin ECEF, originLLA LLA) ECEF {
	return ENUtoECEF(NEDtoENU(p), origin, originLLA)
}
