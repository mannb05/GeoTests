using GeoStats
using GeoIO
import CairoMakie as CMke

image = GeoIO.load("data/sentinel.tif")

boundingbox(image.geometry)
crs(image)

img = image |> Proj(utmsouth(22))
rgb = img |> Select(4 => "R", 3 => "G", 2 => "B")

λ = 5 #intensity param
ascolor(r,g,b) = RGB(λ * r, λ * g, λ * b)

color = rgb |> Map(["R","G","B"] => ascolor => "RGB")

isvegetation(ind) = ind.application_domain == "vegetation"
inds = filter(isvegetation, spectralindices())

n = img |> Select(8 => "N")

ndvi = [rgb n] |> SpectralIndex("NDVI")
mgrvi = [rgb n] |> SpectralIndex("MGRVI")
si = [rgb n] |> SpectralIndex("SI")
spec = [ndvi mgrvi si]

# fig = CMke.Figure()
# ax1 = CMke.Axis(fig[1,1], title="RGB")
# ax2 = CMke.Axis(fig[1,2], title="NDVI")
# ax3 = CMke.Axis(fig[2,1], title="MGRVI")
# ax4 = CMke.Axis(fig[2,2], title="SI")
# viz!(ax1, color.geometry, color = color.RGB)
# viz!(ax2, spec.geometry, color = spec.NDVI)
# viz!(ax3, spec.geometry, color = spec.MGRVI)
# viz!(ax4, spec.geometry, color = spec.SI)
# fig

q30 = quantile(mgrvi.MGRVI, 0.3)
isinside(x) = x < q30
binary = mgrvi |> Map("MGRVI" => isinside => "label")

mask = binary |> ModeFilter()
mask |> viewer

# First objective of assigning true/false value 
# to pixels wrt agriculture achieved

# Now we need to isolate and measure the 
# geometry of the agricultural fields

potrace = mask |> Potrace("label")
region = potrace.geometry[findfirst(potrace.label)]

polys = parent(region)
field = argmax(area, polys)

areafield = area(field) |> u"ha"
perimeterfield = perimeter(field) |> u"km"
println(areafield)
println(perimeterfield)