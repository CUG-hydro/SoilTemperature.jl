"""
计算每层的土壤的热通量 [W m-2]

- `ibeg`: 第一层土壤温度观测具有较大的误差，因此不使用第一层土壤温度
"""
function soil_HeatFlux!(F::V, T::V, κ::V, z::V, z₊ₕ::V;
  F0::FT=NaN, Tsurf::FT=NaN, method="Tsurf", ibeg::Int=1) where {FT<:Real,V<:AbstractVector{FT}}

  N = length(T)
  if method == "Tsurf"
    if ibeg == 1
      F0 = -κ[1] * (Tsurf - T[1]) / (0 - z[1])
    elseif ibeg > 1
      F0 = -κ[ibeg] * (Tsurf - T[ibeg]) / (z[ibeg-1] - z[ibeg])
    end
  end

  @inbounds for i in ibeg:N-1
    d1 = z[i] - z₊ₕ[i]
    d2 = z₊ₕ[i] - z[i+1]
    κ₊ₕ = κ[i] * κ[i+1] * (d1 + d2) / (κ[i] * d2 + κ[i+1] * d1) # Eq. 5.16, 
    # κ₊ₕ = (κ[i] + κ[i+1]) / 2
    Δz₊ₕ = z[i] - z[i+1]
    F[i] = -κ₊ₕ * (T[i] - T[i+1]) / Δz₊ₕ
  end
  F[N] = 0
  return F0
end


"""
包含不参与计算的土壤层

# Arguments
- `method`: 
  + `Tsurf` : Tsurf boundary condition, 第一类边界条件
  + `F0`  : F0 boundary condition, 第二类边界条件
"""
function TsoilEquation(dT, T, soil::Soil, t; method="Tsurf", ibeg::Int=1)
  soil.timestep += 1
  # TODO: 根据t，更新Tsurf
  (; N, Δz, z, z₊ₕ, F, Tsurf) = soil
  (; κ, cv) = soil.param
  F0 = soil_HeatFlux!(F, T, κ, z, z₊ₕ; Tsurf, F0=soil.F0, method, ibeg)

  dT[ibeg] = -(F0 - F[ibeg]) / (Δz[ibeg] * cv[ibeg])
  @inbounds for i in ibeg+1:N
    dT[i] = -(F[i-1] - F[i]) / (Δz[i] * cv[i])
  end
end

# 去掉不参与计算的土壤层
# dT: [ibeg:N]
# T : [ibeg:N]
function TsoilEquation_partial(dT, T, p::Soil, t; method="Tsurf", ibeg::Int=1)
  p.du[ibeg:end] .= dT
  p.u[ibeg:end] .= T
  TsoilEquation(p.du, p.u, p, t; method, ibeg)
  dT .= p.du[ibeg:end]
  return nothing
end

export TsoilEquation, TsoilEquation_partial, soil_HeatFlux!
