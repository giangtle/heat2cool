function [air_1, air_3, air_4, air_5, H2O_consumption] = cal_air_loop_condition(T_env, T_cold_aisle, T_hot_aisle, H_room)
% air_4 stream is conditioned based on the cold aisle temperature and required humidity.
air_4 = struct(); air_4.phase = "gas"; air_4.T = T_cold_aisle;
air_4.yH2O = H_room*cal_yH2Osat(T_cold_aisle);
air_4.yO2 = 0.21*(1-air_4.yH2O);
air_4.yN2 = 0.79*(1-air_4.yH2O);
air_4.n = 1;
% Calculate air_3 streams and water consumption:
[air_3, H2O_consumption] = cal_cooling_tower(air_4, T_env);

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