using GeoStats
using GeoIO
import CairoMakie as CMke

image = GeoIO.load("data/sentinel.tif")

boundingbox(image.geometry)
crs(image)