function [air_out_avg, air_out_dynamic, T_bed_f, t_ads] = cal_adsorption_column(air_in, T_bed_0, mSi, yH2O_ads, capacity)
% Adsorption process.
[~, liq_stream] = cal_isenthalpic_condensation(air_in);
if isstruct(liq_stream)
    fprintf("ERROR: input stream (in) is not all gaseous\n")
    return
end
if ~(isfield(air_in, "yH2O"))
    fprintf("ERROR: input stream (in) does not contain water")
end
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

% Needed input: column Silica mass - m_Si %[g]
mH2O_end = mSi*capacity;
tspan = [0 (mSi*capacity/(air_in.yH2O*air_in.n*H2O.MW*0.01))];
y0 = [0, T_bed_0];
opt = odeset('Events', @myEvent);
[t, y] = ode15s(@(t, y) find_prime(t, y, air_in), tspan, y0, opt);
% adsorption time:
t_ads = t(end);

% Store T_bed_final:
T_bed_f = y(end,2);
% Store exit stream condition:
T_bed_m = y(:,2);
air_out_dynamic = cal_out_stream(air_in, T_bed_m, yH2O_ads);
air_out_dynamic.t = t;

% out average stream:
air_out_avg = cal_out_stream(air_in, T_bed_m, yH2O_ads);
air_out_avg.T = trapz(t, y(:,2))/t_ads;

    function dy = find_prime(~, y, in)
    % Unknowns - y: 1. mH2O, bed     2. T_bed
    dy = zeros(size(y));
    T_bed = y(2);
    out = cal_out_stream(in, T_bed, yH2O_ads);
    dH = cal_stream_enthalpy(only_air(in))-cal_stream_enthalpy(only_air(out));
    % Mass of water:
    dy(1) = ( in.n*in.yH2O - out.n*out.yH2O)*H2O.MW;
    % (mSi*Si.Cp + mH2O*H2O.Cp)*dT/dt = in.yH2O*in.n*H2O.MW*Q_adsorption
    dy(2) = ( dy(1)*Q_adsorption + dH )/( y(1)*H2O.Cp + mSi*Si.Cp );
    end

    function [value, isterminal, direction] = myEvent(~, y)
    value      = y(1) - mH2O_end;
    isterminal = 1;   % Stop the integration
    direction  = 0;
    end
end

function out = cal_out_stream(in, T_bed, yH2O_ads)
% exit air stream:
    out = struct();
    out.phase = "gas";
    out.yH2O = yH2O_ads;
    out.n = in.n*(1-in.yH2O)./(1-out.yH2O);
    % air stream compositions:
    gas_list = {'yCH4', 'yCO', 'yCO2', 'yH2', 'yO2', 'yN2'}; % yH2O not in list.
    for i = 1:length(gas_list)
        gas = gas_list{i};
        if isfield(in, gas)
            out.(gas) = in.(gas)*in.n./out.n;
        end
    end
    out.T = T_bed;
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