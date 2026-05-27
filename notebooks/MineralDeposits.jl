using Revise
using GeoTests

using GeoStats
using GeoIO
using PairPlots
import GLMakie as Mke

url = "https://zenodo.org/record/7051975/files/drillholes.csv?download=1"
csv = download(url, tempname()*".csv")
dtable = GeoIO.load(csv, coords = ("X", "Y", "Z"))

selectholdid = Select("HOLEID")
selectgrades = Select(
    "Cu ppm" => "Cu",
    "Au ppm" => "Au",
    "Ag ppm" => "Ag",
    "S ppm" => "S",
) → Functional(x -> 1e-4*x*u"percent")

dclean = selectholdid ⊔ selectgrades
dtable = dclean(dtable)
# dtable |> Select("Cu","Au","Ag","S") |> values |> pairplot

bbox = boundingbox(dtable.geometry)
bsize = (25.0u"m", 25.0u"m", 12.5u"m")
grid = CartesianGrid(extrema(bbox)..., bsize)

# viz(dtable.geometry, color = "black")
# viz!(grid, alpha = 0.2)
# Mke.current_figure()

shadow(point) = point |> Shadow("xy")
points = shadow.(dtable.geometry)
chull = convexhull(points)

# viz(chull)
# viz!(points, color = "black")
# Mke.current_figure()

active = findall(h -> shadow(centroid(h)) ∈ chull, grid)
blocks = view(grid, active)

zcoord(point) = coords(point).z

ztable = @chain dtable begin
    @groupby(:HOLEID)
    @transform(:z = zcoord(:geometry), :geometry = shadow(:geometry))
    @combine(:z = first(:z), :geometry = first(:geometry))
end

centroids = unique(shadow.(centroid.(blocks)))
ztable = ztable |> Select("z") |> Interpolate(centroids, model=IDW())

p(h) = shadow(centroid(h))
z(h) = zcoord(centroid(h))

zdict = Dict(ztable.geometry .=> ztable.z)
active = findall(h -> z(h) < zdict[p(h)], blocks)
blocks = view(blocks, active)

grades = dtable |> Select("Cu","Au","Ag","S")

preproc = CLR() → ProjectionPursuit()
samples, cache = apply(preproc, grades)

maxlag = 300.0u"m"
vs = setdiff(names(samples), ["geometry"])
gs = [EmpiricalVariogram(samples, v, maxlag = maxlag) for v in vs]
γs = [GeoStatsFunctions.fit(Variogram, g, h -> 1 / h^2) for g in gs]

function gammaplot(n, g, γ)
    fig = Mke.Figure()
    Mke.Axis(fig[1,1], title = n)
    funplot!(fig, g, maxlag = maxlag)
    funplot!(fig, γ, maxlag = maxlag)
end

# gammaplot(vs[1], gs[1], γs[1])

interps = map(vs, γs) do v, γ
    samples |> Select(v) |> InterpolateNeighbors(blocks, model = Kriging(γ))
end

interp = reduce(hcat, interps)

estim = revert(preproc, interp, cache)
# estim |> Select("Ag") |> viewer

μ = mean(estim.Cu) - 0.1
σ = std(estim.Cu)

f(Cu) = 1 / (1 + exp(-(Cu - μ) / σ))

recov = estim |> Map("Cu" => f => "f")

ton = 1000u"kg"
kg = 1000u"g"

ρ = 2.75 * ton / 1u"m^3"
P = 4000 / ton
Cₘ = 4 / ton
Cₚ = 10 / ton

model = @transform([estim recov],
    :value = volume(:geometry) * ρ * ((:Cu / 100) * :f * P - (Cₘ + Cₚ))
)

model |> Select("Cu") |> viewer