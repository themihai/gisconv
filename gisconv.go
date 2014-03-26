package gisconv

import (
	"math"
)

const (
	radToDeg = 180 / math.Pi
	degToRad = math.Pi / 180

	a    = 6377563.396
	b    = 6356256.909  // Airy 1830 major & minor semi-axes
	f0   = 0.9996012717 // NatGrid scale factor on central meridian
	lat0 = 49 * degToRad
	lon0 = -2 * degToRad // NatGrid true origin
	n0   = -100000.0
	e0   = 400000.0        // northing & easting of true origin, metres
	e2   = 1 - (b*b)/(a*a) // eccentricity squared
	n    = (a - b) / (a + b)
	n2   = n * n
	n3   = n * n * n
)

// func OsGridToLatLong takes the northing and easting params and returns latitude and longitude
func OsGridToLatLong(northing, easting float64) (float64, float64) {
	lat := lat0
	m := 0.0
	for northing-n0-m >= 1e-5 { // until < 0.01mm
		lat = (northing-n0-m)/(a*f0) + lat
		ma := (1 + n + (5/4)*n2 + (5/4)*n3) * (lat - lat0)
		mb := (3*n + 3*n*n + (21/8)*n3) * math.Sin(lat-lat0) * math.Cos(lat+lat0)
		mc := ((15/8)*n2 + (15/8)*n3) * math.Sin(2*(lat-lat0)) * math.Cos(2*(lat+lat0))
		md := (35 / 24) * n3 * math.Sin(3*(lat-lat0)) * math.Cos(3*(lat+lat0))
		m = b * f0 * (ma - mb + mc - md) // meridional arc
	}
	cosLat := math.Cos(lat)
	sinLat := math.Sin(lat)
	nu := a * f0 / math.Sqrt(1-e2*sinLat*sinLat)                 // transverse radius of curvature
	rho := a * f0 * (1 - e2) / math.Pow(1-e2*sinLat*sinLat, 1.5) // meridional radius of curvature
	eta2 := nu/rho - 1
	tanLat := math.Tan(lat)
	tan2lat := tanLat * tanLat
	tan4lat := tan2lat * tan2lat
	tan6lat := tan4lat * tan2lat
	secLat := 1 / cosLat
	nu3 := nu * nu * nu
	nu5 := nu3 * nu * nu
	nu7 := nu5 * nu * nu
	vii := tanLat / (2 * rho * nu)
	viii := tanLat / (24 * rho * nu3) * (5 + 3*tan2lat + eta2 - 9*tan2lat*eta2)
	ix := tanLat / (720 * rho * nu5) * (61 + 90*tan2lat + 45*tan4lat)
	x := secLat / nu
	xi := secLat / (6 * nu3) * (nu/rho + 2*tan2lat)
	xii := secLat / (120 * nu5) * (5 + 28*tan2lat + 24*tan4lat)
	xiia := secLat / (5040 * nu7) * (61 + 662*tan2lat + 1320*tan4lat + 720*tan6lat)
	de := easting - e0
	de2 := de * de
	de3 := de2 * de
	de4 := de2 * de2
	de5 := de3 * de2
	de6 := de4 * de2
	de7 := de5 * de2
	lat = lat - vii*de2 + viii*de4 - ix*de6
	lon := lon0 + x*de - xi*de3 + xii*de5 - xiia*de7
	return lat * radToDeg, lon * radToDeg
}
