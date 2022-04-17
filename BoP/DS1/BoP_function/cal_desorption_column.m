function [exhaust_out, steam_avg, steam_dynamic] = cal_desorption_column(exhaust_in, T_bed_1, t_ads, mSi, T_desorp)
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
% Column saturation point:
capacity = 0.3; % [g H2O/g Si]
% Heat of adsorption:
Q_adsorption = 2800;    % [J/g]
% Consistency check:
if (exhaust_in.T <= T_bed_1 || exhaust_in.T <= T_desorp)
    fprintf("ERROR: exhaust temperature too low.")
    return
end

% Initial water mass:
mH2O_0 = mSi*capacity;
% Unknowns: 1: T_out
options = optimset("MaxFunEvals", 10000, "MaxIter", 10000, 'Display', 'off');
% initial guess:
T0 = exhaust_in.T-20;
f = @(T) equation_set(T, t_ads);
[T_out, fval] = fsolve(f, T0, options);
if max(abs(fval))>1e-6
    return
end
[t_desorp, t_heat_up, ~] = cal_desorption_time(exhaust_in, T_out);

% Check for thermal crossing:
if T_out < T_desorp
    fprintf("ERROR: exit exhaust temperature calculated to be lower than desorption temperature.\nNot hot enough exhaust.\n")
    return
end

% return solution:
exhaust_out = exhaust_in;
exhaust_out.T = T_out;

% solve streams:
[steam_avg, steam_dynamic] = cal_steam_stream(exhaust_in, T_out, T_desorp, t_desorp, t_heat_up);

    function F = equation_set(T_out, t_ads)
        F = cal_desorption_time(exhaust_in, T_out) - t_ads;
    end

    function [t_desorp, t_heat_up, t_evap] = cal_desorption_time(exhaust_in, T_out)
        dH = cal_exhaust_Q(exhaust_in, T_out);
        t_heat_up = ( ( mSi*Si.Cp + mH2O_0*H2O.Cp )*(T_desorp-T_bed_1) )/dH;
        t_evap = mH2O_0/( dH/Q_adsorption );
        t_desorp = t_heat_up + t_evap;
    end
end

function dH = cal_exhaust_Q(exhaust_in, T_out)
exhaust_out = exhaust_in;
exhaust_out.T = T_out;
dH = cal_stream_enthalpy(exhaust_in) - cal_stream_enthalpy(exhaust_out);
end

function [steam_avg, steam_dynamic] = cal_steam_stream(exhaust_in, T_out, T_desorp, t_desorp, t_heat_up)
Q_adsorption = 2800;    % [J/g]
H2O_MW = 18.01528;      % [g/mol]
dH = cal_exhaust_Q(exhaust_in, T_out);
t1 = linspace(0, t_heat_up, 20).';
n1 = zeros(size(t1));
t2 = linspace(t_heat_up, t_desorp, 20).';
n2 = ones(size(t2))*dH/(Q_adsorption*H2O_MW);
t = [t1; t2];
n = [n1; n2];
n_avg = trapz(t,n)/t_desorp;

% Assign values:
steam_avg = struct();
steam_avg.phase = "gas";
steam_avg.yH2O = 1;
steam_avg.T = T_desorp;

steam_dynamic = steam_avg;
steam_dynamic.n = n; steam_dynamic.t = t;

steam_avg.n = n_avg;
end