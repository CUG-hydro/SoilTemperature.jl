using Plots, Printf
gr(framestyle=:box)

# inner_optimizer = GradientDescent()
# options = Optim.Options(show_trace=true)

Δz = [2.5, 5, 5, 5, 5, 35, 45, 115, 205] ./ 100
z, z₊ₕ, Δz₊ₕ = soil_depth_init(Δz)
# _z = [4, 12, 14, 16, 20, 24, 28, 32, 36, 42, 50, 52] * 25.4/1000 # inch to mm

function plot_soil(i; ibeg=1)
  i2 = i + ibeg - 1
  title = @sprintf("layer %d: depth = %d cm", i2, -z[i2] * 100)

  time_min, time_max = minimum(t), maximum(t)
  ticks = time_min:Dates.Day(7):time_max
  xticks = ticks, Dates.format.(ticks, "mm-dd")

  plot(; title, xticks)
  plot!(t, yobs[:, i], label="OBS")
  plot!(t, ysim[:, i], label="SIM")
end

function init_soil(; Tsurf=20.0, dt=3600.0, soil_type=1, ibeg=2)
  # Δz = fill(0.025, N)
  # Δz = [2.5, 5, 5, 5, 5, 35, 45, 115, 205] ./ 100
  N = length(Δz)
  z, z₊ₕ, Δz₊ₕ = soil_depth_init(Δz)

  m_sat = θ_S[soil_type] * ρ_wat * Δz # kg/m2
  m_ice = 0 * m_sat
  m_liq = 0.8 * m_sat
  Tsoil = deepcopy(Tsoil0)

  κ, cv = soil_properties_thermal(Δz, Tsoil, m_liq, m_ice;
    soil_type, method="apparent-heat-capacity")
  param = SoilParam{Float64}(; N, κ, cv)
  Soil{Float64}(; N, dt, z, z₊ₕ, Δz, Δz₊ₕ, param, Tsurf, Tsoil, ibeg)
end

function goal(theta; method="ODE")
  soil = init_soil(; soil_type=7)
  solver = Tsit5()
  ysim = model_Tsoil_sim(soil, Tsurf, theta; method, solver)
  obs = yobs[:, 2:end][:]
  sim = ysim[:, 2:end][:]
  # of_MSE(obs, sim)
  gof = GOF(obs, sim)
  -gof.NSE
end
