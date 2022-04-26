function BoP = solve_BoP(air, CH4, H2O, W_e, T_env, T_cold_aisle, T_hot_aisle, H_room)
% FC system additional information:
W_loss = 150;               % [Watts]

% Desorption/Adsorption information:
Desiccant = struct();
Desiccant.id = "Si";
Desiccant.m = 5000; % [g]
Desiccant.capacity = 0.2; % [g H2O per g Desiccant]
Desiccant.Cp = 1.13;   % [J/g/K]
Desiccant.Q_adsorption = 2800;    % [J/g]

% Guess of desorption temperature:
T_desorp = 50+273.15;

% Exhaust stream:
ex_FC = cal_FC_black_box(CH4, H2O, air, W_e, W_loss);

% air-1 conditions
[air_1, air_3, air_4, air_5, H2O_consumption] = cal_air_loop_condition(T_env, T_cold_aisle, T_hot_aisle, H_room);
[air_1, air_2, ex_1, T_desorp, t_cycle] = cal_adsorp_desorp_cycle(ex_FC, air_1, Desiccant, T_desorp, air_3.yH2O);
air_loop = cal_air_loop(air_1, air_2, air_3, air_4, air_5, H2O_consumption);

% water recollection:
water_collection = cal_water_recollection(ex_1, T_env);

% Final assignment:
BoP = assign_stream();

    function BoP = assign_stream()
        BoP = struct();
        % ENERGY:
        LHV = 802340;
        BoP.LHV = (CH4.n*LHV);
        BoP.EE = W_e/(CH4.n*LHV);
        BoP.W_e = W_e;
        BoP.W_server = air_loop.W_server;
        BoP.W_reject_1 = air_loop.W_reject_1;
        BoP.W_reject_2 = air_loop.W_reject_2;
        BoP.W_condenser = water_collection.W_condenser;

        % MATERIAL FLOW:
        % MAIN SYSTEM:
        BoP.n_CH4_FC = CH4.n;
        BoP.n_H2O_FC = H2O.n;
        BoP.n_air_FC = air.n;
        BoP.yH2O_air_FC = air.yH2O;
        BoP.yO2_air_FC = air.yO2;
        BoP.yN2_air_FC = air.yN2;
        BoP.n_ex = ex_1.n;
        BoP.yH2O_ex = ex_1.yH2O;
        BoP.yCO2_ex = ex_1.yCO2;
        BoP.yO2_ex = ex_1.yO2;
        BoP.yN2_ex = ex_1.yN2;
        % AIR LOOP:
        BoP.n_air_1 = air_loop.air_1.n;
        BoP.yH2O_air_1 = air_loop.air_1.yH2O;
        BoP.yO2_air_1 = air_loop.air_1.yO2;
        BoP.yN2_air_1 = air_loop.air_1.yN2;
        BoP.H_room = air_loop.air_4.yH2O/(cal_yH2Osat(air_loop.air_4.T));
        % WATER RECOLLECTION:
        BoP.water_consumption_from_air_loop = air_loop.H2O_consumption.n;
        BoP.water_consumption_from_FC = H2O.n;
        BoP.water_collected = water_collection.ex_cond.n;
        BoP.net_water = water_collection.ex_cond.n - air_loop.H2O_consumption.n - H2O.n;

        % TEMPERATURE (in oC):
        % CH4 in SYSTEM:
        BoP.T_CH4_FC = CH4.T - 273.15;
        % H2O in SYSTEM:
        BoP.T_H2O_FC = H2O.T - 273.15;
        % AIR in SYSTEM:
        BoP.T_air_FC = air.T - 273.15;
        % EXHAUST in SYSTEM:
        BoP.T_ex_FC = ex_FC.T - 273.15;
        BoP.T_ex_1 = ex_1.T - 273.15;
        BoP.T_ex_2 = water_collection.ex_2.T - 273.15;
        % AIRLOOP:
        BoP.T_desorp = T_desorp -273.15;
        BoP.T_air_1 = air_loop.air_1.T - 273.15;
        BoP.T_air_2 = air_loop.air_2.T - 273.15;
        BoP.T_air_3 = air_loop.air_3.T - 273.15;
        BoP.T_air_4 = air_loop.air_4.T - 273.15;
        BoP.T_air_5 = air_loop.air_5.T - 273.15;
        BoP.T_H2O_consumption = air_loop.H2O_consumption.T - 273.15;
        BoP.H_air_1 = air_loop.air_1.yH2O/(cal_yH2Osat(air_loop.air_1.T));
        BoP.H_air_2 = air_loop.air_2.yH2O/(cal_yH2Osat(air_loop.air_2.T));
        BoP.H_air_3 = air_loop.air_3.yH2O/(cal_yH2Osat(air_loop.air_3.T));
        BoP.H_air_4 = air_loop.air_4.yH2O/(cal_yH2Osat(air_loop.air_4.T));
        BoP.H_air_5 = air_loop.air_5.yH2O/(cal_yH2Osat(air_loop.air_5.T));
        
        % Adsorption - desorption cycle:
        BoP.mDesiccant = Desiccant.m;
        BoP.capacity = Desiccant.capacity;
        BoP.t_cycle = t_cycle;
    end
end

function water_collection = cal_water_recollection(ex_1, T_env)
ex_2 = ex_1; ex_2.T = T_env + 5;
[ex_2, ex_cond] = cal_isothermal_condensation(ex_2);
W_condenser = cal_stream_enthalpy(ex_1) - cal_stream_enthalpy(ex_2) - cal_stream_enthalpy(ex_cond);
% Report back:
water_collection = struct();
water_collection.ex_2 = ex_2;
water_collection.ex_cond = ex_cond;
water_collection.W_condenser = W_condenser;
end

function air_loop = cal_air_loop(air_1, air_2, air_3, air_4, air_5, H2O_consumption)
air_3.n = air_2.n;
air_4.n = air_1.n;
air_5.n = air_1.n;
H2O_consumption.n = air_4.n-air_3.n;
% Calculate energy streams W_server and W_air_HX:
W_server = - ( cal_stream_enthalpy(air_4) - cal_stream_enthalpy(air_5) );
W_reject_1 = cal_stream_enthalpy(air_2) - cal_stream_enthalpy(air_3);
W_reject_2 = cal_stream_enthalpy(air_5) - cal_stream_enthalpy(air_1);
% save info to air loop struct:
air_loop = struct();
air_loop.W_server = W_server; air_loop.W_reject_1 = W_reject_1; air_loop.W_reject_2 = W_reject_2;
air_loop.T_cold_aisle = air_4.T;
air_loop.T_hot_aisle = air_5.T;
air_loop.T_H2O = H2O_consumption.T;
air_loop.air_1 = air_1;
air_loop.air_2 = air_2;
air_loop.air_3 = air_3;
air_loop.air_4 = air_4;
air_loop.air_5 = air_5;
air_loop.H2O_consumption = H2O_consumption;
end