package geodesy

import "math"

// XTE2D computes signed cross-track error from point P to infinite line A->B.
// Uses only horizontal (E/N or X/Y) components, ignoring altitude.
//
// Sign convention:
//   - positive: P is left of track A->B
//   - negative: P is right of track A->B
//   - zero: P is on the track
//
// Returns distance in same units as input (meters for ENU).
func XTE2D(A, B, P ENU) float64 {
	ABx := B.X - A.X
	ABy := B.Y - A.Y

	APx := P.X - A.X
	APy := P.Y - A.Y

	return (ABx*APy - ABy*APx) / math.Hypot(ABx, ABy)
}

// XTE2DSegment computes signed cross-track error from point P to segment A->B.
// If projection falls outside the segment, returns distance to nearest endpoint.
//
// Use this when trajectory may overshoot segment endpoints.
func XTE2DSegment(A, B, P ENU) float64 {
	ABx := B.X - A.X
	ABy := B.Y - A.Y

	APx := P.X - A.X
	APy := P.Y - A.Y

	ABLen2 := ABx*ABx + ABy*ABy
	if ABLen2 == 0 {
		return math.Hypot(APx, APy)
	}

	// Projection factor
	t := (APx*ABx + APy*ABy) / ABLen2

	switch {
	case t < 0:
		return math.Hypot(APx, APy)
	case t > 1:
		return math.Hypot(P.X-B.X, P.Y-B.Y)
	default:
		return (ABx*APy - ABy*APx) / math.Sqrt(ABLen2)
	}
}

// XTE3D computes unsigned cross-track error from point P to infinite line A->B in 3D.
// Uses all three components (E/N/U or X/Y/Z).
//
// Note: Result is always positive (no left/right sign in 3D).
//
// Returns distance in same units as input (meters for ENU).
func XTE3D(A, B, P ENU) float64 {
	AB := Vec3{X: B.X - A.X, Y: B.Y - A.Y, Z: B.Z - A.Z}
	AP := Vec3{X: P.X - A.X, Y: P.Y - A.Y, Z: P.Z - A.Z}

	C := cross(AP, AB)
	return norm(C) / norm(AB)
}

// XTE3DSegment computes unsigned minimum distance from point P to segment A->B in 3D.
// If projection falls outside the segment, returns distance to nearest endpoint.
func XTE3DSegment(A, B, P ENU) float64 {
	AB := Vec3{X: B.X - A.X, Y: B.Y - A.Y, Z: B.Z - A.Z}
	AP := Vec3{X: P.X - A.X, Y: P.Y - A.Y, Z: P.Z - A.Z}

	ABLen2 := AB.X*AB.X + AB.Y*AB.Y + AB.Z*AB.Z
	if ABLen2 == 0 {
		return norm(AP)
	}

	// Projection factor
	t := (AP.X*AB.X + AP.Y*AB.Y + AP.Z*AB.Z) / ABLen2

	switch {
	case t < 0:
		return norm(AP)
	case t > 1:
		BP := Vec3{X: P.X - B.X, Y: P.Y - B.Y, Z: P.Z - B.Z}
		return norm(BP)
	default:
		C := cross(AP, AB)
		return norm(C) / math.Sqrt(ABLen2)
	}
}

// cross computes 3D cross product of vectors A and B.
func cross(A, B Vec3) Vec3 {
	return Vec3{
		X: A.Y*B.Z - A.Z*B.Y,
		Y: A.Z*B.X - A.X*B.Z,
		Z: A.X*B.Y - A.Y*B.X,
	}
}

// norm computes magnitude of vector V.
func norm(V Vec3) float64 {
	return math.Sqrt(V.X*V.X + V.Y*V.Y + V.Z*V.Z)
}
