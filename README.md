# Soil Temperature

<!-- [![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://jl-pkgs.github.io/SoilTemperature.jl/stable) -->
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://jl-pkgs.github.io/SoilTemperature.jl/dev)
[![CI](https://github.com/jl-pkgs/SoilTemperature.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/jl-pkgs/SoilTemperature.jl/actions/workflows/CI.yml)
[![codecov](https://codecov.io/gh/jl-pkgs/SoilTemperature.jl/branch/master/graph/badge.svg?token=gNwrAxE8oz)](https://codecov.io/gh/jl-pkgs/SoilTemperature.jl/tree/master/src)

> 求解土壤热通量方程

<!-- - [x] `Bonan 2021`：计算速度快，但公式复杂
- [x] `diffeq`: 公式清晰，但计算速度过慢 -->

![](image/case01_Tsoil_CUG_ODE.png)
图1. 大气系土壤温度观测，以及土壤温度模拟结果。

<!-- ## 求解方案 -->
<!-- > 注意，sink需要划分到每一层的蒸发量 -->

## References

- <https://github.com/jl-pkgs/HydroTools.jl/blob/master/src/Soil/soil_moisture.jl>

- <https://github.com/amireson/RichardsEquation/blob/master/Richards%20Equation.ipynb>


<!-- ## 测试站点

<https://mesonet.agron.iastate.edu/agclimate/hist/hourly.php> -->
