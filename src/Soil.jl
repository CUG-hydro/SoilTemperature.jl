export SoilParam, Soil
using Parameters
using Printf

export Soil, SoilParam, ParamVanGenuchten

## 结构体形式的参数
abstract type AbstractSoilParam{FT} end

@with_kw mutable struct ParamVanGenuchten{T} <: AbstractSoilParam{T}
  θ_sat::T = 0.287       # [m3 m-3]
  θ_res::T = 0.075       # [m3 m-3]
  Ksat::T = 34 / 3600    # [cm s-1]
  α::T = 0.027
  n::T = 3.96
  m::T = 1.0 - 1.0 / n
end

# 参数优化过程中，可能需要优化的参数
# 一个重要的经验教训，不要去优化`m`，NSE会下降0.2
@with_kw mutable struct SoilParam{FT}
  ## Parameter: 土壤水力
  N::Int = 10
  same_layer = false

  ## Parameter: 土壤热力
  κ::Vector{FT} = fill(2.0, N)         # thermal conductivity [W m-1 K-1]
  cv::Vector{FT} = fill(2.0 * 1e6, N)  # volumetric heat capacity [J m-3 K-1]
end


@with_kw_noshow mutable struct Soil{FT}
  N::Int = 10                        # layers of soil
  ibeg::Int = 1                      # index of the first layer，边界层条件指定
  inds_obs::Vector{Int} = ibeg:N     # indices of observed layers

  dt::Float64 = 3600                 # 时间步长, seconds
  z::Vector{FT} = zeros(FT, N)       # m, 向下为负
  z₊ₕ::Vector{FT} = zeros(FT, N)
  Δz::Vector{FT} = zeros(FT, N)
  Δz₊ₕ::Vector{FT} = zeros(FT, N)

  # 温度
  Tsoil::Vector{FT} = fill(NaN, N)   # [°C]
  κ₊ₕ::Vector{FT} = zeros(FT, N - 1)  # thermal conductivity at interface [W m-1 K-1]
  F::Vector{FT} = zeros(FT, N)       # heat flux, [W m-2]
  Tsurf::FT = FT(NaN)                  # surface temperature, [°C]
  F0::FT = FT(NaN)                   # heat flux at the surface, [W m-2]，向下为负
  G::FT = FT(NaN)                    # [W m-2]，土壤热通量

  ## Parameter: [水力] + [热力]参数
  param::SoilParam{FT} = SoilParam{FT}(; N)
  # param_water::ParamVanGenuchten{FT} = ParamVanGenuchten{FT}()

  # ODE求解临时变量
  u::Vector{FT} = fill(NaN, N)  # [°C], 为了从ibeg求解地温，定义的临时变量
  du::Vector{FT} = fill(NaN, N) # [°C]

  timestep::Int = 0                  # 迭代次数
end

function Soil(Δz::Vector{FT}; kw...) where {FT}
  N = length(Δz)
  z, z₊ₕ, Δz₊ₕ = soil_depth_init(Δz)
  soil = Soil{Float64}(; N, z, z₊ₕ, Δz, Δz₊ₕ, kw...)
  # update K and ψ
  cal_K!(soil)
  cal_ψ!(soil)
  return soil
end

# θ = fill(0.1, N)
# ψ = van_Genuchten_ψ.(θ; param=param_water)
# θ0 = 0.267
# ψ0 = van_Genuchten_ψ(θ0; param=param_water)
# dt = 5 # [s]
# sink = ones(N) * 0.3 / 86400 # [cm s⁻¹], 蒸发速率

function Base.show(io::IO, param::SoilParam{T}) where {T<:Real}
  # (; use_m, same_layer) = param
  printstyled(io, "Parameters: \n", color=:blue, bold=true)
  # println("[use_m = $use_m, same_layer = $same_layer]")
  println(io, "-----------------------------")
  print_var(io, param, :κ)
  print_var(io, param, :cv; scale=1e6)
  println(io, "-----------------------------")
  return nothing
end


function Base.show(io::IO, x::Soil{T}) where {T<:Real}
  param = x.param

  printstyled(io, "Soil{$T}: ", color=:blue)
  printstyled(io, "N = $(x.N), ibeg=$(x.ibeg), ", color=:blue, underline=true)
  print_index(io, x.inds_obs; prefix="inds_obs =")

  printstyled(io, "Soil Temperature: \n", color=:blue, bold=true)
  print_var(io, x, :Tsoil)
  print_var(io, x, :Tsurf)
  
  show(io, param)
  return nothing
end

function print_selected(io::IO, name::String, method::String)
  if name[1:5] == method[1:5]
    printstyled(io, "   [$name]\n", bold=true, color=:green)
  else
    printstyled(io, "   [$name]\n", bold=true)
  end
end

function print_var(io::IO, x, var; scale=nothing, digits=3, color=:blue, used=true)
  value = getfield(x, var)
  name = @sprintf("%-5s", string(var))
  _color = used ? color : :white
  printstyled(io, " - $name: "; color=_color)
  if isnothing(scale)
    println(io, round.(value; digits))
  else
    println(io, "$(round.(value/scale; digits)) * $scale")
  end
end

function print_index(io::IO, inds; prefix="", color=:blue, underline=true)
  if length(unique(diff(inds))) == 1
    n = length(inds)
    printstyled(io, "$prefix $(inds[1]):$(inds[end]) [n=$n] \n"; color, underline)
  else
    printstyled(io, "$prefix $inds \n"; color, underline)
  end
end
