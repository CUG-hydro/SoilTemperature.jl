using SoilTemperature, Plots, Test, RTableTools, Dates
using OrdinaryDiffEq
import HydroTools: sceua, GOF, of_KGE, of_NSE
using Optim
includet("main_Tsoil.jl")

d = fread("$dir_soil/data/CUG_TS_202306.csv")
t = d.time
A = Matrix(d[:, 2:end]) #|> drop_missing 

ibeg = 2
Tsurf = A[:, ibeg]
yobs = A[:, ibeg:end]
Tsoil0 = A[1, :]

solver = Tsit5()
# solver = Rosenbrock23()
# solver = Rodas5(autodiff=false)

begin
  soil = init_soil(; soil_type=7)
  method = "Bonan"
  # method = "ODE"

  if method == "Bonan"
    ysim = solve_Tsoil_Bonan(soil, Tsurf)
  elseif method == "ODE"
    ysim = solve_Tsoil_ODE(soil, Tsurf; solver)
  end

  plot(
    [plot_soil(i; ibeg) for i in 1:8]...,
    size=(1200, 800),
  )
  # savefig("case01_Tsoil_CUG_$(method).png")
end

begin
  soil = init_soil(; soil_type=7)
  x0 = [soil.κ; soil.cv]
  lower = [fill(0.1, 9); fill(0.1, 9) * 1e6]
  upper = [fill(10.0, 9); fill(5.0, 9) * 1e6]

  f(theta) = goal(theta)
  @time theta, feval, exitflag = sceua(f, x0, lower, upper; maxn=Int(5 * 1e4))
end

ysim = model_Tsoil_sim(soil, Tsurf, theta;)

# r = optimize(f, lower, upper, x0, Fminbox(inner_optimizer), options)
# theta = r.minimizer
