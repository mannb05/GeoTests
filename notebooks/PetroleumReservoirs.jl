using GeoStats
using GeoIO
import GLMakie as Mke

norne₁ = GeoIO.load("data/norne1.vtu")
norne₂ = GeoIO.load("data/norne2.vtu")

clean = Select(
  "porosity" => "ϕ",
  "saturation_oil" => "So",
  "saturation_gas" => "Sg",
  "saturation_water" => "Sw",
  "density_oil" => "ρo [kg/m^3]",
  "density_gas" => "ρg [kg/m^3]",
  "density_water" => "ρw [kg/m^3]"
) → Unitify()

reservoir₁ = clean(norne₁)
reservoir₂ = clean(norne₂)

mass(reservoir) = @transform(reservoir,
    :MOIP = :ρo * :So * :ϕ * volume(:geometry),
    :MGIP = :ρg * :Sg * :ϕ * volume(:geometry),
    :MWIP = :ρw * :Sw * :ϕ * volume(:geometry)
)

mass₁ = mass(reservoir₁)
mass₂ = mass(reservoir₂)

zones = @transform(mass₁, :MHIP = :MOIP + :MGIP) |>
        Select("MHIP") |> GHC(3, 500u"m", nmax=1000) |> Rename("label" => "zone")

carbon₁ = mass₁ |> Select("MOIP" => "MOIP₁", "MGIP" => "MGIP₁")
carbon₂ = mass₂ |> Select("MOIP" => "MOIP₂", "MGIP" => "MGIP₂")

data = [carbon₁ carbon₂ zones]

summary = @chain data begin
    @groupby(:zone)
    @transform(:delta = :MOIP₁ + :MGIP₁ - :MOIP₂ - :MGIP₂)
    @combine(:depletion = sum(:delta))
end

summary |> Unit("depletion" => u"Mg")