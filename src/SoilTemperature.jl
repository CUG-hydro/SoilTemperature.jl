module SoilTemperature

using DiffEqBase
import HydroTools: sceua, GOF, of_KGE, of_NSE
using Parameters

# greet() = print("Hello World!")
include("Soil.jl")
include("EquationTsoil.jl")
include("Solve_Tsoil.jl")
include("Soil_depth.jl")
include("Soil_properties_thermal.jl")

dir_soil = "$(@__DIR__)/.." |> abspath
export dir_soil

# export soil_temperature!, soil_temperature_F0!
export TsoilEquation, TsoilEquation_partial

end # module SoilTemperature
