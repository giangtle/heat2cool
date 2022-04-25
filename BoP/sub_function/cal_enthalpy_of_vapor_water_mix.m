function dH = cal_enthalpy_of_vapor_water_mix(in)
if in.phase ~= "liq-vap"
    fprintf("ERROR: this function is only to calculate stream enthalpy of mixed saturated vapor and liquid water stream.");
    return
end
h_sat_vap = enthalpy("H2O", 100 + 273.15);
h_sat_liq = enthalpy("H2Ol", 100 + 273.15);
dH = in.n*( in.q*h_sat_vap + (1-in.q)*h_sat_liq );
end