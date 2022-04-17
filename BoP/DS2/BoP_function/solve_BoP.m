function BoP = solve_BoP(air, CH4, n_H2O, W_e, T_env, T_cold_aisle, T_hot_aisle, yH2O_ads, capacity)
% FC system additional information:
W_loss = 150;               % [Watts]
% Desorption/Adsorption column additional information:
mSi = 1000;                % [g]
% Initial guess of desorption temperature:
T_desorp = 40+273.15;

% Water supply stream:
H2O = struct(); H2O.phase = "liq"; H2O.yH2Ol = 1; H2O.T = T_env; H2O.n = n_H2O;

% exhaust calculation:
ex_FC = cal_FC_black_box(CH4, H2O, air, W_e, W_loss);

% air-1 conditions
air_1 = cal_air_loop_condition(T_env, T_cold_aisle, T_hot_aisle, yH2O_ads);
[air_1, air_2, ex_1, T_desorp, t_cycle] = cal_ads_desorption_cycle(ex_FC, air_1, mSi, T_desorp, yH2O_ads, capacity);
air_loop = cal_air_loop(air_1, air_2, T_env, T_cold_aisle, T_hot_aisle, yH2O_ads);

% water recollection:
water_collection = cal_water_recollection(ex_1, T_env);

% Final assignment:
BoP = assign_stream();

    function BoP = assign_stream()
        BoP = struct();
        % ENERGY:
        LHV = 802340; EE = W_e/(CH4.n*LHV);
        BoP.EE = EE;
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
        BoP.n_air_loop_H2O_consumption = air_loop.H2O_consumption.n;
        BoP.H_room = air_loop.air_4.yH2O/(cal_water_y_sat(air_loop.air_4.T));
        % WATER RECOLLECTION:
        BoP.n_flue_condensate = water_collection.n_flue_condensate;
        BoP.net_water = BoP.n_flue_condensate - BoP.n_air_loop_H2O_consumption;

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
        BoP.H_air_1 = air_loop.air_1.yH2O/(cal_water_y_sat(air_loop.air_1.T));
        BoP.H_air_2 = air_loop.air_2.yH2O/(cal_water_y_sat(air_loop.air_2.T));
        BoP.H_air_3 = air_loop.air_3.yH2O/(cal_water_y_sat(air_loop.air_3.T));
        BoP.H_air_4 = air_loop.air_4.yH2O/(cal_water_y_sat(air_loop.air_4.T));
        BoP.H_air_5 = air_loop.air_5.yH2O/(cal_water_y_sat(air_loop.air_5.T));
        
        % Adsorption - desorption cycle:
        BoP.mSi = mSi;
        BoP.capacity = capacity;
        BoP.t_cycle = t_cycle;
    end
end

function water_collection = cal_water_recollection(ex_1, T_env)
ex_2 = ex_1; ex_2.T = T_env + 5;
[gas_stream, liq_stream] = cal_isothermal_condensation(ex_2);
W_condenser = cal_stream_enthalpy(ex_1) - cal_stream_enthalpy(gas_stream) - cal_stream_enthalpy(liq_stream);
n_flue_condensate = liq_stream.n;
% Report back:
water_collection = struct();
water_collection.ex_2 = ex_2;
water_collection.gas_stream = gas_stream;
water_collection.liq_stream = liq_stream;
water_collection.W_condenser = W_condenser;
water_collection.n_flue_condensate = n_flue_condensate;
end

function air_loop = cal_air_loop(air_1, air_2, T_env, T_cold_aisle, T_hot_aisle, yH2O_ads)
[~, air_3, air_4, air_5, H2O_consumption] = cal_air_loop_condition(T_env, T_cold_aisle, T_hot_aisle, yH2O_ads);
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
air_loop.T_cold_aisle = T_cold_aisle;
air_loop.T_hot_aisle = T_hot_aisle;
air_loop.T_H2O = T_env;
air_loop.air_1 = air_1;
air_loop.air_2 = air_2;
air_loop.air_3 = air_3;
air_loop.air_4 = air_4;
air_loop.air_5 = air_5;
air_loop.H2O_consumption = H2O_consumption;
end