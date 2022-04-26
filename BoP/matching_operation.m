function BoP = matching_operation(T_env, H_env, T_cold_aisle, T_hot_aisle, H_room, CH4, H2O, air)
%if stream variables don't exist use the base case conditions: 
if exist("CH4", "var") && exist("H2O", "var") && exist("air", "var")
    fprintf("Using user-decalred operation.\n");
else
    fprintf("Using base case operation.\n");
    [CH4, H2O, air] = base_case(T_env, H_env);
end

% Initial guess:
We_0 = 1000;

% Solve:
f = @(W_e) eqn_matching_We_and_Wserver(air, CH4, H2O, W_e, T_env, T_cold_aisle, T_hot_aisle, H_room);
options = optimset("MaxFunEvals", 10000, "MaxIter", 10000, 'Display', 'off');
[W_e, fval] = fsolve(f, We_0, options);
if max(abs(fval))>1e-6
    return
end
BoP = solve_BoP(air, CH4, H2O, W_e, T_env, T_cold_aisle, T_hot_aisle, H_room);
end

function F = eqn_matching_We_and_Wserver(air, CH4, H2O, W_e, T_env, T_cold_aisle, T_hot_aisle, H_room)
BoP = solve_BoP(air, CH4, H2O, W_e, T_env, T_cold_aisle, T_hot_aisle, H_room);
F = BoP.W_e - BoP.W_server;
end