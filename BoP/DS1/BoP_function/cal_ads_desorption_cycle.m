function [steam_prod, air_1, air_2, ex_1] = cal_ads_desorption_cycle(ex_FC, air_1, T_bed_0, mSi, T_desorp, dT_out, yH2O_ads)
% Initial guess:
y0 = 0.4;
% Condition:
T_ex_out = T_desorp + dT_out;

% Solve what molar flow rate need to be that ex_2.T meets condition:
options = optimset("MaxFunEvals", 10000, "MaxIter", 10000, 'Display', 'off');
f = @(y) equation_set(y, air_1, ex_FC, T_bed_0, mSi, T_ex_out, T_desorp, yH2O_ads);
[n_air_1, fval] = fsolve(f, y0, options);

if max(abs(fval))>1e-6
    return
end

[~, air_1, air_2, steam_prod, ex_1] = cal_desorption_exhaust_temp(n_air_1, air_1, ex_FC, T_bed_0, mSi, T_desorp, yH2O_ads);

end

function F = equation_set(y, air_in, ex_in, T_bed_0, mSi, T_ex_out, T_desorp, yH2O_ads)
F = cal_desorption_exhaust_temp(y, air_in, ex_in, T_bed_0, mSi, T_desorp, yH2O_ads) - T_ex_out;
end

function [T, air_1, air_2, steam_prod, ex_1] = cal_desorption_exhaust_temp(y, air_1, ex_FC, T_bed_0, mSi, T_desorp, yH2O_ads)
air_1.n = y;
[air_2, ~, T_bed_f, t_ads] = cal_adsorption_column(air_1, T_bed_0, mSi, yH2O_ads);
[ex_1, steam_prod, ~] = cal_desorption_column(ex_FC, T_bed_f, t_ads, mSi, T_desorp);
T = ex_1.T;
end