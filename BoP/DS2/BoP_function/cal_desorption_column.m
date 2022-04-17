function [ex_out, T_desorp, t_desorp] = cal_desorption_column(ex_in, T_bed_i, mSi, T_desorp, capacity)
% Constraints/Assumptions:
%   _ integral rate of water removal from t=0 to t=t_ads equals total water
%   (mH2O_0)
%   _ Fixed exit exhaust temperature.
% Unknowns:
%   _ exhaust_out.T

% Known constants:
% Adsorbent: Silica
Si = struct();
Si.Cp = 1.13;   % [J/g/K]
% Water properties:
H2O = struct();
H2O.Cp = 4.18;  % [J/g/K]
H2O.MW = 18.01528;  % [g/mol]
% Heat of adsorption:
Q_adsorption = 2800;    % [J/g]

% Water mole fraction of exhaust out is silica gel equilibrium at T_desorp:
% If exhaust out has higher water mole fraction, T_desorp needs to increase
[T_desorp, ex_out] = cal_T_desorp(ex_in, T_desorp);

% Consistency check:
if (ex_in.T <= T_bed_i || ex_in.T <= T_desorp)
    fprintf("ERROR: exhaust temperature too low.")
    return
end

% Initial water mass:
mH2O_0 = mSi*capacity;

% dH:
dH_heat_up = cal_dH_heat_up(ex_in, T_desorp);
dH_evap = cal_stream_enthalpy(only_air(ex_in))-cal_stream_enthalpy(only_air(ex_out));
% heat up time:
t_heat_up = ( ( mSi*Si.Cp + mH2O_0*H2O.Cp )*(T_desorp-T_bed_i) )/dH_heat_up;

% evaporation time:
t_evap = mH2O_0*Q_adsorption/dH_evap;

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