function [ex_1, T_desorp, t_desorp] = cal_desorption_column(ex_FC, T_bed_i, Desiccant, T_desorp)
% Constraints/Assumptions:
%   _ integral rate of water removal from t=0 to t=t_ads equals total water
%   (mH2O_0)
%   _ Fixed exit exhaust temperature.
% Unknowns:
%   _ exhaust exit temperature: ex_1.T

% Water properties:
H2O = struct();
H2O.Cp = 4.18;  % [J/g/K]
H2O.MW = 18.01528;  % [g/mol]

% Water mole fraction of exhaust out is silica gel equilibrium at T_desorp:
% If exhaust out has higher water mole fraction, T_desorp needs to increase
[T_desorp, ex_1] = cal_T_desorp(ex_FC, T_desorp);

% Consistency check:
if (ex_FC.T <= T_bed_i || ex_FC.T <= T_desorp)
    fprintf("ERROR: exhaust temperature too low.")
    return
end

% Initial water mass:
mH2O_0 = Desiccant.m*Desiccant.capacity;

% dH:
dH_heat_up = cal_dH_heat_up(ex_FC, T_desorp);
dH_evap = cal_stream_enthalpy(only_air(ex_FC))-cal_stream_enthalpy(only_air(ex_1));

% dQ:
Q_heat_up = ( ( Desiccant.m*Desiccant.Cp + mH2O_0*H2O.Cp )*(T_desorp-T_bed_i) );
Q_desorp = mH2O_0*Desiccant.Q_adsorption;

% heat up time:
t_heat_up = Q_heat_up/dH_heat_up;

% evaporation time:
t_evap = Q_desorp/dH_evap;

% Total desorption time:
t_desorp = t_heat_up + t_evap;
end

function out = only_air(in)
% air only:
    out = in;
    out.n = in.n*(1-in.yH2O);
    % air part compositions:
    gas_list = {'yCH4', 'yCO', 'yCO2', 'yH2', 'yO2', 'yN2'}; % yH2O not in list.
    for i = 1:length(gas_list)
        gas = gas_list{i};
        if isfield(in, gas)
            out.(gas) = in.(gas)*in.n/out.n;
        end
    end
    out.yH2O = 0;
end

function dH_heat_up = cal_dH_heat_up(ex_in, T_desorp)
ex_out = ex_in;
ex_out.T = T_desorp;
dH_heat_up = cal_stream_enthalpy(ex_in) - cal_stream_enthalpy(ex_out);
end