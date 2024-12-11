export θ_S, ρ_wat


const ρ_wat = 1000.0                       # Density of water, [kg/m3]
const ρ_ice = 917.0                        # Density of ice, [kg/m3]
const λ_fus = 0.3337 * 1e6                 # Heat of fusion for water at 0℃, [J/kg]
const K0 = 273.15                          # zero degree Celsius, [K]
const tfrz = K0                            # freezing Temperature (K)

## Those values are not used in MOST time
"""
# Initialize soil texture variables

## Soil texture classes (Cosby et al. 1984)
%  1: sand
%  2: loamy sand
%  3: sandy loam
%  4: silty loam
%  5: loam
%  6: sandy clay loam
%  7  silty clay loam
%  8: clay loam
%  9: sandy clay
% 10: silty clay
% 11: clay

## References
- Cosby et al. 1984. Water Resources Research 20:682-690, Soil texture classes
- Clapp and Hornberger. 1978. Water Resources Research 14:601-604
"""
#  (Cosby et al. 1984. Water Resources Research 20:682-690)
SILT = [5.0, 12.0, 32.0, 70.0, 39.0, 15.0, 56.0, 34.0, 6.0, 47.0, 20.0] # Percent silt
SAND = [92.0, 82.0, 58.0, 17.0, 43.0, 58.0, 10.0, 32.0, 52.0, 6.0, 22.0] # Percent sand
CLAY = [3.0, 6.0, 10.0, 13.0, 18.0, 27.0, 34.0, 34.0, 42.0, 47.0, 58.0] # Percent clay

# Volumetric soil water content (%) at saturation (porosity)
# (Clapp and Hornberger. 1978. Water Resources Research 14:601-604)
θ_S = [0.395, 0.410, 0.435, 0.485, 0.451, 0.420, 0.477, 0.476, 0.426, 0.492, 0.482]


"""
    soil_properties_thermal(dz::AbstractVector, Tsoil::AbstractVector,
        m_liq::AbstractVector, m_ice::AbstractVector;
        soil_type::Integer=1, method="excess-heat")

# Arguments

- `dz`    : the thickness of each soil layer (m)
- `m_liq` : Unfrozen water, liquid (kg H2O/m2)
- `m_ice` : Frozen water, ice (kg H2O/m2)
- `Tsoil` : Soil temperature of each soil layer (℃)

- `method`: method of phase change
- `soil_type`: 
  + `1`: sand

# Return

- `κ` : thermal conductivity, [W/m/K]
- `cv`: heat capacity, [J/m3/K]
"""
function soil_properties_thermal(dz::AbstractVector, Tsoil::AbstractVector,
  m_liq::AbstractVector, m_ice::AbstractVector;
  soil_type::Integer=1, method="excess-heat")

  # Volumetric soil water content (%) at saturation (porosity)
  # (Clapp and Hornberger. 1978. Water Resources Research 14:601-604)
  θ_S = [0.395, 0.410, 0.435, 0.485, 0.451, 0.420, 0.477, 0.476, 0.426, 0.492, 0.482]
  
  _c_wat = 4188.0                         # Specific heat of water (J/kg/K)
  _c_ice = 2117.27                        # Specific heat of ice (J/kg/K)

  cv_wat = _c_wat * ρ_wat # Heat capacity of water (J/m3/K)
  cv_ice = _c_ice * ρ_ice # Heat capacity of ice (J/m3/K)

  tfrz = 0.0                            # Freezing point of water [k]

  ## --- Physical constants in physcon structure
  κ_wat = 0.57                          # Thermal conductivity of water (W/m/K)
  κ_ice = 2.29                          # Thermal conductivity of ice (W/m/K)

  ## --- Model run control parameters
  n = length(dz)
  κ = zeros(n)
  cv = zeros(n)

  k = soil_type
  @fastmath @inbounds for i = 1:n
    # --- Volumetric soil water and ice
    θ_liq = m_liq[i] / (ρ_wat * dz[i])
    θ_ice = m_ice[i] / (ρ_ice * dz[i])

    # Fraction of total volume that is liquid water
    fᵤ = θ_liq / (θ_liq + θ_ice)

    # --- Dry thermal conductivity (W/m/K) from 
    ρ_b = 2700 * (1 - θ_S[k]) # density of soil soliads, bulk density (kg/m3)
    κ_dry = (0.135 * ρ_b + 64.7) / (2700 - 0.947 * ρ_b)  # Eq. 5.27

    # --- Kersten number and unfrozen and frozen values
    S_e = min((θ_liq + θ_ice) / θ_S[k], 1) # Soil water relative to saturation
    Ke_f = S_e

    if (SAND[k] < 50)
      Ke_u = 1 + log10(max(S_e, 0.1))
    else
      Ke_u = 1 + 0.7 * log10(max(S_e, 0.05))
    end
    Ke = Tsoil[i] >= tfrz ? Ke_u : Ke_f

    ## --- Soil solids thermal conducitivty (W/m/K)
    q = SAND[k] / 100         # Quartz fraction
    κ_o = q > 0.2 ? 2.0 : 3.0 # Thermal conductivity of other minerals (W/m/K)
    κ_q = 7.7                 # Thermal conductivity of q (W/m/K)

    # Thermal conductivity of soil solids (W/m/K)
    κ_sol = κ_q^q * κ_o^(1 - q)  # Eq. 5.31

    # --- Saturated thermal conductivity (W/m/K) and unfrozen and frozen values
    κ_sat = κ_sol^(1 - θ_S[k]) * κ_wat^(fᵤ * θ_S[k]) * κ_ice^((1 - fᵤ) * θ_S[k]) # Eq. 5.30

    κ_sat_u = κ_sol^(1 - θ_S[k]) * κ_wat^θ_S[k] # Eq. 5.28
    κ_sat_f = κ_sol^(1 - θ_S[k]) * κ_ice^θ_S[k] # Eq. 5.29

    # --- Thermal conductivity (W/m/K) and unfrozen and frozen values
    κ[i] = (κ_sat - κ_dry) * Ke + κ_dry

    κ_u = (κ_sat_u - κ_dry) * Ke_u + κ_dry
    κ_f = (κ_sat_f - κ_dry) * Ke_f + κ_dry

    ## --- Heat capacity of soil solids (J/m3/K)
    cv_sol = 1.926e6

    # --- Heat capacity (J/m3/K) and unfrozen and frozen values
    cv[i] = (1 - θ_S[k]) * cv_sol + cv_wat * θ_liq + cv_ice * θ_ice # Eq. 5.32

    cv_u = (1 - θ_S[k]) * cv_sol + cv_wat * (θ_liq + θ_ice)
    cv_f = (1 - θ_S[k]) * cv_sol + cv_ice * (θ_liq + θ_ice)

    # --- Adjust heat capacity and thermal conductivity if using apparent heat capacity
    if method == "apparent-heat-capacity"
      # 这里考虑了结冰和融化的过程
      tinc = 0.5 # Temperature range for freezing and thawing [k]
      # Heat of fusion (J/m3), equivalent to 
      # ql = λ_fus * (m_liq + m_ice) / dz
      ql = λ_fus * (ρ_wat * θ_liq + ρ_ice * θ_ice)

      # Heat capacity and thermal conductivity, Eq. 5.39
      if Tsoil[i] > tfrz + tinc
        cv[i] = cv_u
        κ[i] = κ_u
      elseif tfrz - tinc <= Tsoil[i] <= tfrz + tinc
        cv[i] = (cv_f + cv_u) / 2 + ql / (2 * tinc)
        κ[i] = κ_f + (κ_u - κ_f) * (Tsoil[i] - tfrz + tinc) / (2 * tinc)
      elseif Tsoil[i] < tfrz - tinc
        cv[i] = cv_f
        κ[i] = κ_f
      end
    end
  end
  κ, cv
end

export soil_properties_thermal
