function BoP = solve_BoP(air, CH4, n_H2O, W_e, T_env, T_cold_aisle, T_hot_aisle, T_desorp, yH2O_ads)
% FC system additional information:
W_loss = 150;               % [Watts]
% Desorption/Adsorption column additional information:
mSi = 1000;                % [g]
T_bed_0 = T_desorp;

% Air loop conditions
air_1 = cal_air_loop_condition(T_env, T_cold_aisle, T_hot_aisle, yH2O_ads);
% Water supply stream:
H2O_supply = struct(); H2O_supply.phase = "liq"; H2O_supply.yH2Ol = 1; H2O_supply.T = T_env;
% Surrogate flow of H2O to calculate ex-1 molar information:
H2O = H2O_supply; H2O.n = n_H2O;
[ex_1, ~] = cal_FC_black_box(CH4, H2O, air, W_e, W_loss);
ex_1 = rmfield(ex_1, "T");
% Solve the system through Newton loop:
    T_guess = 400 + 273.15;
    f = cal_ex_FC_T(T_guess, air, CH4, T_env, n_H2O, ex_1, air_1, T_bed_0, mSi, W_e, W_loss, T_desorp, T_cold_aisle, T_hot_aisle, yH2O_ads);
    iteration = 1;
    while abs(f)>1e-3 && iteration < 15
        T_guess = T_guess - f;
        f = cal_ex_FC_T(T_guess, air, CH4, T_env, n_H2O, ex_1, air_1, T_bed_0, mSi, W_e, W_loss, T_desorp, T_cold_aisle, T_hot_aisle, yH2O_ads);
    end
    [~, HXN, air_loop, steam_prod, H2O_FC, ex_FC, ex_1, water_collection] = cal_ex_FC_T(T_guess, air, CH4, T_env, n_H2O, ex_1, air_1, T_bed_0, mSi, W_e, W_loss, T_desorp, T_cold_aisle, T_hot_aisle, yH2O_ads);
% Final assignment:
BoP = assign_stream(air, CH4, H2O_FC, steam_prod, ex_FC, ex_1, HXN, air_loop, W_e, water_collection);
end

function BoP = assign_stream(air, CH4, H2O_FC, steam_prod, ex_FC, ex_1, HXN, air_loop, W_e, water_collection)
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
BoP.n_H2O_FC = H2O_FC.n;
BoP.n_air_FC = air.n;
BoP.yH2O_air_FC = air.yH2O;
BoP.yO2_air_FC = air.yO2;
BoP.yN2_air_FC = air.yN2;
BoP.n_ex = ex_1.n;
BoP.yH2O_ex = ex_1.yH2O;
BoP.yCO2_ex = ex_1.yCO2;
BoP.yO2_ex = ex_1.yO2;
BoP.yN2_ex = ex_1.yN2;
BoP.n_ex_HX1 = HXN.HX1.ex_HX1_in.n;
BoP.n_ex_HX2 = HXN.HX2.ex_HX2_in.n;
BoP.n_steam_prod = steam_prod.n;
if steam_prod.n < H2O_FC.n
    BoP.n_H2O_supplement = H2O_FC.n - steam_prod.n;
    BoP.n_H2O_HX0 = "NA";
elseif HXN.config == "2"
    BoP.n_H2O_supplement = "NA";
    BoP.n_H2O_HX0 = HXN.HX0.H2O_HX0_in.n;
    if HXN.HX0.H2O_HX0_out.phase == "liq-vap"
        BoP.q_H2O_HX0_out = HXN.HX0.H2O_HX0_out.q;
    else
        BoP.q_H2O_HX0_out = "NA";
    end
end
% AIR LOOP:
BoP.n_air_1 = air_loop.air_1.n;
BoP.yH2O_air_1 = air_loop.air_1.yH2O;
BoP.yO2_air_1 = air_loop.air_1.yO2;
BoP.yN2_air_1 = air_loop.air_1.yN2;
BoP.n_air_loop_H2O_consumption = air_loop.H2O_consumption.n;
BoP.H_room = air_loop.air_4.yH2O/(cal_water_y_sat(air_loop.air_4.T));
% WATER RECOLLECTION:
BoP.n_steam_condensate = water_collection.n_steam_condensate;
BoP.n_flue_condensate = water_collection.n_flue_condensate;
BoP.n_water_recollected = water_collection.n_net;
BoP.net_water = water_collection.n_steam_condensate + water_collection.n_flue_condensate - BoP.n_air_loop_H2O_consumption;

% TEMPERATURE (in oC):
% CH4 in SYSTEM:
BoP.T_CH4 = CH4.T - 273.15;
BoP.T_CH4_FC = HXN.HX2.CH4_FC.T - 273.15;
% H2O in SYSTEM:
BoP.T_H2O_FC = H2O_FC.T - 273.15;
BoP.T_steam_prod = steam_prod.T - 273.15;
if HXN.config == "2"
    BoP.T_H2O_HX0_out = HXN.HX0.H2O_HX0_out.T - 273.15;
else
    BoP.T_H2O_HX0_out = "NA";
end
% AIR in SYSTEM:
BoP.T_air = air.T - 273.15;
if HXN.config == "2"
    BoP.T_air_HX0_out = HXN.HX0.air_HX0_out.T - 273.15;
else
    BoP.T_air_HX0_out = "NA";
end
BoP.T_air_FC = HXN.HX1.air_FC.T - 273.15;
% EXHAUST in SYSTEM:
BoP.T_ex_FC = ex_FC.T - 273.15;
BoP.T_ex_1 = ex_1.T - 273.15;
BoP.T_ex_HX1_out = HXN.HX1.ex_HX1_out.T - 273.15;
BoP.T_ex_HX2_out = HXN.HX2.ex_HX2_out.T - 273.15;
% AIRLOOP:
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
end

function [error, HXN, air_loop, steam_prod, H2O_FC, actual_ex_FC, ex_1, water_collection] = cal_ex_FC_T(T_guess, air, CH4, T_env, n_H2O, ex_1, air_1, T_bed_0, mSi, W_e, W_loss, T_desorp, T_cold_aisle, T_hot_aisle, yH2O_ads)
dT_out = 10;
ex_FC = ex_1;
ex_FC.T = T_guess;
[steam_prod, air_1, air_2, ex_1] = cal_ads_desorption_cycle(ex_FC, air_1, T_bed_0, mSi, T_desorp, dT_out, yH2O_ads);
air_loop = cal_air_loop(air_1, air_2, T_env, T_cold_aisle, T_hot_aisle, yH2O_ads);
% Determine what H2O stream should be.
if steam_prod.n < n_H2O
    n_H2O_supply = n_H2O - steam_prod.n;
    H2O_FC = cal_stream_from_mixing_water_steam(n_H2O_supply, T_env, steam_prod);
    [HXN, CH4_FC, air_FC] = case_1_HXN(air, CH4, ex_1);
else
    H2O_FC = steam_prod;
    H2O_FC.n = n_H2O;
    H2O_HX0_in = steam_prod; H2O_HX0_in.n = steam_prod.n - n_H2O;
    [HXN, CH4_FC, air_FC] = case_2_HXN(air, CH4, ex_1, H2O_HX0_in);
end
% Condenser calculation to recollect water:
water_collection = cal_water_recollection(HXN, T_env);
% Difference between guess of T_ex and actual T_ex:
actual_ex_FC = cal_FC_black_box(CH4_FC, H2O_FC, air_FC, W_e, W_loss);
T_ex = actual_ex_FC.T;
error = ex_FC.T - T_ex;
end

function water_collection = cal_water_recollection(HXN, T_env)
W_condenser = 0;
n_steam_condensate = 0;
n_flue_condensate = 0;
if HXN.config == "2"
    % Energy to condense water from HX0:
    water_in = HXN.HX0.H2O_HX0_out;
    if water_in.phase == "liq-vap"
        water_out = water_in; water_out = rmfield(water_out, "q"); water_out.phase = "liq"; water_out.yH2Ol = 1; water_out.T = T_env + 5;
        W_condenser = W_condenser + cal_enthalpy_of_vapor_water_mix(water_in) - cal_stream_enthalpy(water_out);
    else
        water_out = water_in; water_out.phase = "liq"; water_out.yH2Ol = 1; water_out.T = T_env + 5;
        W_condenser = W_condenser + cal_stream_enthalpy(water_in) - cal_stream_enthalpy(water_out);
    end
    n_steam_condensate = n_steam_condensate + water_in.n;
end
% Water collection from HX1 and HX2:
    % HX1:
    in = HXN.HX1.ex_HX1_out;
    out = in; out.T = T_env + 5;
    [gas_stream, liq_stream] = cal_isothermal_condensation(out);
    W_condenser = W_condenser + ( cal_stream_enthalpy(in)-cal_stream_enthalpy(gas_stream)-cal_stream_enthalpy(liq_stream) );
    n_flue_condensate = n_flue_condensate + liq_stream.n;
    % HX2:
    in = HXN.HX2.ex_HX2_out;
    out = in; out.T = T_env + 5;
    [gas_stream, liq_stream] = cal_isothermal_condensation(out);
    W_condenser = W_condenser + ( cal_stream_enthalpy(in)-cal_stream_enthalpy(gas_stream)-cal_stream_enthalpy(liq_stream) );
    n_flue_condensate = n_flue_condensate + liq_stream.n;
% Report back:
water_collection = struct();
water_collection.W_condenser = W_condenser;
water_collection.n_steam_condensate = n_steam_condensate;
water_collection.n_flue_condensate = n_flue_condensate;
water_collection.n_net = n_steam_condensate + n_flue_condensate;
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

function [HXN, CH4_FC, air_FC] = case_1_HXN(air, CH4, ex_1)
% HX1 calculations:
    % Temperature approach:
    dT_ap_HX1 = 20;
    air_FC = air; air_FC.T = ex_1.T - dT_ap_HX1;
    ex_HX1_in = ex_1; 
    ex_HX1_out = ex_HX1_in; ex_HX1_out.T = air.T + dT_ap_HX1;
    % Check condensation:
    ex_HX1_out = check_condensation(ex_HX1_out);
    % Flow rate:
    dQ_HX1_cold = cal_stream_enthalpy(air) - cal_stream_enthalpy(air_FC);
    dQ_HX1_hot = cal_stream_enthalpy(ex_HX1_in) - cal_stream_enthalpy(ex_HX1_out);
    ex_HX1_in.n = ex_1.n*( -dQ_HX1_cold/dQ_HX1_hot );
    ex_HX1_out.n = ex_HX1_in.n;
% HX2 calculations:
    % Flow rate:
    ex_HX2_in = ex_1; ex_HX2_in.n = ex_1.n-ex_HX1_in.n;
    % Temperature approach:
    dT_ap_HX2 = 10;
    [ex_HX2_out, CH4_FC, dT_ap_HX2] = cal_HX_streams(ex_HX2_in, CH4, dT_ap_HX2);
    % Check condensation:
    [~, liq_stream] = cal_isenthalpic_condensation(ex_HX2_out);
    if isstruct(liq_stream)
        if liq_stream.n > 0.001*ex_HX2_out.n
            fprintf("ERROR: exhaust condenses running out of the HX2.\n")
            return
        end
    end
% Final HXN results:
HXN = struct();
HXN.config = "1";
HXN.HX1 = struct();
HXN.HX1.dT_ap = dT_ap_HX1;
HXN.HX1.air = air; HXN.HX1.air_FC = air_FC; HXN.HX1.ex_HX1_in = ex_HX1_in; HXN.HX1.ex_HX1_out = ex_HX1_out;
HXN.HX2 = struct();
HXN.HX2.dT_ap = dT_ap_HX2;
HXN.HX2.CH4 = CH4; HXN.HX2.CH4_FC = CH4_FC; HXN.HX2.ex_HX2_in = ex_HX2_in; HXN.HX2.ex_HX2_out = ex_HX2_out;
end

function [HXN, CH4_FC, air_FC] = case_2_HXN(air, CH4, ex_1, H2O_HX0_in)
flag = 0;
% HX0 calculations:
    % Temperature approach:
    % Assumption: steam condensate at 100oC.
    H2O_HX0_out = struct();
    H2O_HX0_out.phase = "liq";
    H2O_HX0_out.n = H2O_HX0_in.n;
    H2O_HX0_out.yH2Ol = 1;
    H2O_HX0_out.T = 100 + 273.15;
    dQ_HX0_hot = cal_stream_enthalpy(H2O_HX0_in) - cal_stream_enthalpy(H2O_HX0_out);
    % Find hot exit temperature - air_HX0_out:
    H_air_HX0_out = cal_stream_enthalpy(air) + dQ_HX0_hot;
    air_HX0_out = cal_T_from_enthalpy(air, H_air_HX0_out);
    % Check temperature restriction:
    if air_HX0_out.T > (H2O_HX0_in.T-10)
        clear H2O_HX0_out air_HX0_out
        air_HX0_out = air; air_HX0_out.T = H2O_HX0_in.T - 10;
        H2O_HX0_out = struct();
        H2O_HX0_out.phase = "liq-vap";
        H2O_HX0_out.n = H2O_HX0_in.n;
        H2O_HX0_out.T = 100 + 273.15;
        % calculate the stream quality:
        h_H2O_HX0_out = ( cal_stream_enthalpy(H2O_HX0_in) + cal_stream_enthalpy(air) - cal_stream_enthalpy(air_HX0_out) )/( H2O_HX0_in.n );
        h_sat_vap = enthalpy("H2O", 100 + 273.15);
        h_sat_liq = enthalpy("H2Ol", 100 + 273.15);
        H2O_HX0_out.q = ( h_H2O_HX0_out-h_sat_liq )/( h_sat_vap-h_sat_liq );
        flag = 1;
    end
% HX1 calculations:
    % Temperature approach:
    dT_ap_HX1 = 10;
    air_FC = air_HX0_out; air_FC.T = ex_1.T - dT_ap_HX1;
    ex_HX1_in = ex_1; 
    ex_HX1_out = ex_HX1_in; ex_HX1_out.T = air_HX0_out.T + dT_ap_HX1;
    % Check condensation:
    ex_HX1_out = check_condensation(ex_HX1_out);
    % Flow rate:
    dQ_HX1_cold = cal_stream_enthalpy(air_HX0_out) - cal_stream_enthalpy(air_FC);
    dQ_HX1_hot = cal_stream_enthalpy(ex_HX1_in) - cal_stream_enthalpy(ex_HX1_out);
    ex_HX1_in.n = ex_1.n*( -dQ_HX1_cold/dQ_HX1_hot );
    ex_HX1_out.n = ex_HX1_in.n;
% HX2 calculations:
    % Flow rate:
    ex_HX2_in = ex_1; ex_HX2_in.n = ex_1.n-ex_HX1_in.n;
    % Temperature approach:
    dT_ap_HX2 = 10;
    [ex_HX2_out, CH4_FC, dT_ap_HX2] = cal_HX_streams(ex_HX2_in, CH4, dT_ap_HX2);
    % Check condensation:
    [~, liq_stream] = cal_isenthalpic_condensation(ex_HX2_out);
    if isstruct(liq_stream)
        if liq_stream.n > 0.001*ex_HX2_out.n
            fprintf("ERROR: exhaust condenses running out of the HX2.\n")
            return
        end
    end
% Final HXN results:
HXN = struct();
HXN.config = "2";
HXN.HX0 = struct();
if flag == 0
    HXN.HX0.dT_ap = cal_HX_temp_ap(H2O_HX0_in, H2O_HX0_out, air, air_HX0_out);
end
HXN.HX0.air = air; HXN.HX0.air_HX0_out = air_HX0_out; HXN.HX0.H2O_HX0_in = H2O_HX0_in; HXN.HX0.H2O_HX0_out = H2O_HX0_out;
HXN.HX1 = struct();
HXN.HX1.dT_ap = dT_ap_HX1;
HXN.HX1.air_HX0_out = air_HX0_out; HXN.HX1.air_FC = air_FC; HXN.HX1.ex_HX1_in = ex_HX1_in; HXN.HX1.ex_HX1_out = ex_HX1_out;
HXN.HX2 = struct();
HXN.HX2.dT_ap = dT_ap_HX2;
HXN.HX2.CH4 = CH4; HXN.HX2.CH4_FC = CH4_FC; HXN.HX2.ex_HX2_in = ex_HX2_in; HXN.HX2.ex_HX2_out = ex_HX2_out;
end

function stream = check_condensation(stream)
% Check condensation:
[~, liq_stream] = cal_isenthalpic_condensation(stream);
if isstruct(liq_stream)
    if liq_stream.n > 0.001*stream.n
        fprintf("ERROR: exhaust condenses running out of the HX1.\n");
        stream.T = stream.T + 5;
        stream = check_condensation(stream);
    end
end
end