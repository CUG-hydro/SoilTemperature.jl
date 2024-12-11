export model_Tsoil_sim, Tsoil_theta2param, of_MSE
export solve_Tsoil_ODE, solve_Tsoil_Bonan

function Tsoil_theta2param(theta)
  N = length(theta) ÷ 2
  κ = theta[1:N]
  cv = theta[N+1:end]
  return κ, cv
end

function model_Tsoil_sim(soil, Tsurf, theta; method="Bonan", kw...)
  κ, cv = Tsoil_theta2param(theta)
  soil.param.κ .= κ
  soil.param.cv .= cv

  if method == "Bonan"
    ysim = solve_Tsoil_Bonan(soil, Tsurf;)
  elseif method == "ODE"
    ysim = solve_Tsoil_ODE(soil, Tsurf; kw...)
  end
  ysim
end


function solve_Tsoil_Bonan(soil::Soil{FT}, Tsurf::AbstractVector{FT}) where {FT<:Real}
  (; N, inds_obs, ibeg) = soil

  ntime = length(Tsurf)
  R = zeros(ntime, soil.N - ibeg + 1)
  R[1, :] .= soil.Tsoil[ibeg:end]

  for i = 2:ntime
    soil_temperature!(soil, Tsurf[i]; ibeg)
    R[i, :] .= soil.Tsoil[ibeg:end]
  end

  inds = inds_obs .- ibeg .+ 1
  R[:, inds]
end


"""
    solve_Tsoil_ODE(soil, Tsurf; solver)

solver = Tsit5()
solver = Rosenbrock23()
solver = Rodas5(autodiff=false)  
"""
function solve_Tsoil_ODE(soil, Tsurf; solver, reltol=1e-3, abstol=1e-3, verbose=false)
  (; N, inds_obs, ibeg, dt) = soil

  ntime = length(Tsurf)
  u0 = soil.Tsoil[ibeg:end]

  _TsoilEquation(dT, T, p, t) = TsoilEquation_partial(dT, T, p, t; ibeg)
  tspan = (0, dt)
  prob = ODEProblem(_TsoilEquation, u0, tspan, soil)

  R = zeros(ntime, N - ibeg + 1)
  R[1, :] .= soil.Tsoil[ibeg:end]

  for i = 2:ntime
    soil.Tsurf = Tsurf[i]
    prob.u0 .= soil.Tsoil[ibeg:end]

    sol = solve(prob, solver; reltol, abstol, saveat=dt)
    soil.Tsoil[ibeg:end] .= sol.u[end] # 更新这个时刻的结果
    R[i, :] .= soil.Tsoil[ibeg:end]
  end

  inds = @. inds_obs - ibeg + 1
  R[:, inds]
end
