# geodesy

Go package for geodetic coordinate conversions, distance calculations, and cross-track error (XTE) computation on the WGS-84 ellipsoid.

## Install

```bash
go get github.com/pryamcem/geodesy
```

## Features

- **Coordinate conversions**: LLA (lat/lon/alt) <-> ECEF(Earth-Centered, Earth-Fixed) <-> ENU(East-North-Up)
- **Distance**: Haversine (great-circle) and Vincenty (ellipsoidal)
- **Cross-track error**: 2D/3D XTE to infinite lines and segments

## Basic usage

```go
import "github.com/pryamcem/geodesy"

// Distance between two points
a := geodesy.LLA{Lat: 55.7558, Lon: 37.6173, Alt: 0}
b := geodesy.LLA{Lat: 59.9343, Lon: 30.3351, Alt: 0}

dist := geodesy.VincentyDistance(a, b)   // accurate, ellipsoidal
dist = geodesy.GreatCircleDistance(a, b) // fast, spherical

// Convert to local ENU coordinates
origin := a
p := geodesy.LLAtoENU(b, origin)
// p.X = east (m), p.Y = north (m), p.Z = up (m)

// Cross-track error
A := geodesy.ENU{X: 0, Y: 0, Z: 0}
B := geodesy.ENU{X: 100, Y: 0, Z: 0}
P := geodesy.ENU{X: 50, Y: 5, Z: 0}

xte := geodesy.XTE2D(A, B, P) // signed, horizontal
```
