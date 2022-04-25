function BoP = matching_operation(T_env, H_env, T_cold_aisle, T_hot_aisle, yH2O_ads)
% Flows:
air = struct(); air.phase = "gas"; air.n = SLM_to_mol_sec(130); air.T = T_env;
air.yH2O = H_env*cal_yH2Osat(T_env); air.yO2 = 0.21*(1-air.yH2O); air.yN2 = 0.79*(1-air.yH2O);
CH4 = struct(); CH4.phase = "gas"; CH4.n = SLM_to_mol_sec(3.65); CH4.yCH4 = 1; CH4.T = T_env;
n_H2O = SLM_to_mol_sec(7.3);
% Initial guess:
We_0 = 1000;
% Solve:
f = @(W_e) eqn_matching_We_and_Wserver(air, CH4, n_H2O, W_e, T_env, T_cold_aisle, T_hot_aisle, yH2O_ads);
options = optimset("MaxFunEvals", 10000, "MaxIter", 10000, 'Display', 'off');
[W_e, fval] = fsolve(f, We_0, options);
if max(abs(fval))>1e-6
    return
end
BoP = solve_BoP(air, CH4, n_H2O, W_e, T_env, T_cold_aisle, T_hot_aisle, yH2O_ads);
end

function F = eqn_matching_We_and_Wserver(air, CH4, n_H2O, W_e, T_env, T_cold_aisle, T_hot_aisle, yH2O_ads)
BoP = solve_BoP(air, CH4, n_H2O, W_e, T_env, T_cold_aisle, T_hot_aisle, yH2O_ads);
F = BoP.W_e - BoP.W_server;
end