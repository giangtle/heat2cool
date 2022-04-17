function mixed_stream = cal_stream_from_mixing_water_steam(n_H2O_supply, T_env, steam_prod)
n_total = steam_prod.n + n_H2O_supply;
% Water supply stream:
H2O_supply = struct(); H2O_supply.phase = "liq"; H2O_supply.yH2Ol = 1; H2O_supply.T = T_env; H2O_supply.n = n_H2O_supply;
% Initial guess:
T0 = steam_prod.T;
% Solve:
options = optimset("MaxFunEvals", 10000, "MaxIter", 10000, 'Display', 'on');
f = @(T) eqn_mix_stream_temp(T, n_total, H2O_supply, steam_prod);
[T, fval] = fsolve(f, T0, options);
if max(abs(fval))>1e-6
    return
end
mixed_stream = stream_from_T(T, n_total);
if mixed_stream.T < 100+273.15
    mixed_stream.phase = "liq";
    mixed_stream = rmfield(mixed_stream, "yH2O");
    mixed_stream.yH2Ol = 1;
    mixed_stream.T = T_env;
end
end

function F = eqn_mix_stream_temp(T, n_total, H2O_supply, steam_prod)
mixed_stream = stream_from_T(T, n_total);
F = cal_stream_enthalpy(mixed_stream) - ( cal_stream_enthalpy(H2O_supply)+cal_stream_enthalpy(steam_prod) );
end

function mixed_stream = stream_from_T(T, n_total)
mixed_stream = struct();
mixed_stream.n = n_total;
    mixed_stream.phase = "gas";
    mixed_stream.yH2O = 1;
    mixed_stream.T = T;
end