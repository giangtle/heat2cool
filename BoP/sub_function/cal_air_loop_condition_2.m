function [air_1, air_3, air_4, air_5, H2O_consumption] = cal_air_loop_condition_2(T_env, T_cold_aisle, T_hot_aisle, H_room)
% Temperature air_3 stream is limitted by environment temperature.
air_4 = struct(); air_4.phase = "gas"; air_4.T = T_cold_aisle;
air_4.yH2O = H_room*cal_yH2Osat(T_cold_aisle);
air_4.yO2 = 0.21*(1-air_4.yH2O);
air_4.yN2 = 0.79*(1-air_4.yH2O);
air_4.n = 1;
% Calculate air_3 streams and water consumption:
[H2O_consumption, air_4] = cal_cooling_tower(air_3, T_env, T_cold_aisle);
% Check to see if air_4 is above saturated water mole fraction:
if air_4.yH2O > cal_yH2O_sat(air_4.T)
    fprintf("ERROR: check cal_air_loop_condition.m air_4 yH2O mole fraction exceeds saturated water mole fraction.");
    return
end

% Set stream 5 and 1:
air_1 = air_4; air_5 = air_4;
air_5.T = T_hot_aisle;
% If T_hot_aisle > T_env +5, hot air is cooled through another air heat exchanger:
if T_hot_aisle > T_env+5
    air_1.T = T_env+5;
else
    air_1.T = T_hot_aisle;
end

% Remove all molar flowrate information:
air_1 = rmfield(air_1, "n");
air_3 = rmfield(air_3, "n");
air_4 = rmfield(air_4, "n");
air_5 = rmfield(air_5, "n");
H2O_consumption = rmfield(H2O_consumption, "n");
end

function [H2O_consumption, air_4] = cal_cooling_tower(air_3, T_env, T_cold_aisle)
% Temperature air_4 stream is defined by cold aisle temperature.
air_4 = struct(); air_4.phase = "gas"; air_4.T = T_cold_aisle;
% Water supply stream:
H2O_consumption = struct(); H2O_consumption.phase = "liq"; H2O_consumption.yH2Ol = 1; H2O_consumption.T = T_env;
% Calculate the mole of water consumption:
% Initial guess:
nH2O_0 = 0;
% Solve:
f = @(nH2O) eqn_cooling_tower(nH2O, H2O_consumption, air_3, air_4);
options = optimset("MaxFunEvals", 10000, "MaxIter", 10000, 'Display', 'off');
[nH2O, fval] = fsolve(f, nH2O_0, options);
if max(abs(fval))>1e-6
    return
end
[H2O_consumption, air_4] = cal_streams(nH2O, H2O_consumption, air_3, air_4);
end

function F = eqn_cooling_tower(nH2O, H2O_supply, air_3, air_4)
[H2O_supply, air_4] = cal_streams(nH2O, H2O_supply, air_3, air_4);
F = cal_stream_enthalpy(air_3) + cal_stream_enthalpy(H2O_supply) - cal_stream_enthalpy(air_4);
end

function [H2O_supply, air_4] = cal_streams(nH2O, H2O_supply, air_3, air_4)
% air mass balance:
air_4.n = air_3.n + nH2O;
air_4.yH2O = (nH2O + air_3.n*air_3.yH2O)/air_4.n;
air_4.yO2 = air_3.n*air_3.yO2/air_4.n;
air_4.yN2 = air_3.n*air_3.yN2/air_4.n;
% water supply mole:
H2O_supply.n = nH2O;
end