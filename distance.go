package geodesy

import "math"

// GreatCircleDistance calculates the shortest distance between two points on
// a sphere using the Haversine formula. Uses mean Earth radius.
//
// Returns distance in meters.
//
// Note: Faster but less accurate than VincentyDistance for long distances.
// https://en.wikipedia.org/wiki/Great-circle_distance
func GreatCircleDistance(p1, p2 LLA) float64 {
	lat1 := deg2rad(p1.Lat)
	lon1 := deg2rad(p1.Lon)
	lat2 := deg2rad(p2.Lat)
	lon2 := deg2rad(p2.Lon)

	dLat := lat2 - lat1
	dLon := lon2 - lon1

	sinDLat := math.Sin(dLat / 2)
	sinDLon := math.Sin(dLon / 2)

	a := sinDLat*sinDLat +
		math.Cos(lat1)*math.Cos(lat2)*sinDLon*sinDLon

	c := 2 * math.Atan2(math.Sqrt(a), math.Sqrt(1-a))

	return EarthRadius * c
}

// VincentyDistance calculates the geodesic distance between two points on the
// WGS-84 ellipsoid using Vincenty's iterative formula.
//
// Returns distance in meters.
//
// Note: More accurate than GreatCircleDistance, especially for long distances,
// as it accounts for Earth's ellipsoidal shape.
// https://en.wikipedia.org/wiki/Vincenty's_formulae
func VincentyDistance(p1, p2 LLA) float64 {
	lat1 := deg2rad(p1.Lat)
	lon1 := deg2rad(p1.Lon)
	lat2 := deg2rad(p2.Lat)
	lon2 := deg2rad(p2.Lon)

	const maxIter = 200
	const tol = 1e-12

	U1 := math.Atan((1 - wgs84F) * math.Tan(lat1))
	U2 := math.Atan((1 - wgs84F) * math.Tan(lat2))
	L := lon2 - lon1

	lambda := L
	var sinSigma, cosSigma, sigma float64
	var sinAlpha, cos2Alpha, cos2SigmaM float64

	for range maxIter {
		sinLambda := math.Sin(lambda)
		cosLambda := math.Cos(lambda)

		sinSigma = math.Sqrt(
			math.Pow(math.Cos(U2)*sinLambda, 2) +
				math.Pow(
					math.Cos(U1)*math.Sin(U2)-
						math.Sin(U1)*math.Cos(U2)*cosLambda,
					2),
		)

		if sinSigma == 0 {
			return 0 // coincident points
		}

		cosSigma = math.Sin(U1)*math.Sin(U2) +
			math.Cos(U1)*math.Cos(U2)*cosLambda

		sigma = math.Atan2(sinSigma, cosSigma)

		sinAlpha = math.Cos(U1) * math.Cos(U2) * sinLambda / sinSigma
		cos2Alpha = 1 - sinAlpha*sinAlpha

		if cos2Alpha != 0 {
			cos2SigmaM = cosSigma - 2*math.Sin(U1)*math.Sin(U2)/cos2Alpha
		} else {
			cos2SigmaM = 0 // equatorial line
		}

		C := wgs84F / 16 * cos2Alpha * (4 + wgs84F*(4-3*cos2Alpha))

		lambdaPrev := lambda
		lambda = L + (1-C)*wgs84F*sinAlpha*
			(sigma+C*sinSigma*(cos2SigmaM+
				C*cosSigma*(-1+2*cos2SigmaM*cos2SigmaM)))

		if math.Abs(lambda-lambdaPrev) < tol {
			break
		}
	}

	u2 := cos2Alpha * (wgs84A*wgs84A - wgs84B*wgs84B) / (wgs84B * wgs84B)

	A := 1 + u2/16384*(4096+u2*(-768+u2*(320-175*u2)))
	B := u2 / 1024 * (256 + u2*(-128+u2*(74-47*u2)))

	deltaSigma := B * sinSigma *
		(cos2SigmaM + B/4*(cosSigma*(-1+2*cos2SigmaM*cos2SigmaM)-
			B/6*cos2SigmaM*(-3+4*sinSigma*sinSigma)*
				(-3+4*cos2SigmaM*cos2SigmaM)))

	return wgs84B * A * (sigma - deltaSigma)
}
