function [gas_stream, liq_stream] = cal_isenthalpic_condensation(in)
if in.phase == "liq"
    gas_stream = 0;
    liq_stream = in;
else
    yH2O_sat = cal_water_y_sat(in.T);
    if in.yH2O <= yH2O_sat
        gas_stream = in;
        liq_stream = 0;
    else
        [gas_stream, liq_stream] = adiabatic_condesation(in);
    end
end
end

function [gas_stream, liq_stream] = adiabatic_condesation(in)
options = optimset("MaxFunEvals", 10000, "MaxIter", 10000, 'Display', 'off');

% initial guess:
y0 = [0, in.T];

% unknown: n_H2Ol, T
f = @(y) equation_set(y, in);
[y, fval] = fsolve(f,y0,options);

if max(abs(fval))>1e-6
    return
end
[gas_stream, liq_stream] = stream_from_solution(y, in);
end

function F = equation_set(y, in)
[gas_stream, liq_stream] = stream_from_solution(y, in);

% Equations:
% 1: yH2O = yH2O_sat = exp(-gibbs("evap", y(2))/(8.314472*y(2))
F(1) = ( in.n*in.yH2O - y(1) )/( in.n - y(1) ) - cal_water_y_sat(y(2));
% 2: enthalpy(gas_stream) + enthalpy(liq_stream) - enthalpy(in)
F(2) = cal_stream_enthalpy(liq_stream) + cal_stream_enthalpy(gas_stream) - cal_stream_enthalpy(in);
end

function [gas_stream, liq_stream] = stream_from_solution(y, in)
% unknown: n_H2Ol, T
% start out two streams:
liq_stream = struct();
gas_stream = struct();
% mass balance to give composition & # of molar flow in gas stream and liq stream
    % liq_stream:
    liq_stream.n = y(1);
    liq_stream.yH2Ol = 1;
    %gas stream:
    gas_stream.phase = "gas";
    gas_stream.n = in.n - y(1);
    gas_list = {'yCH4', 'yCO', 'yCO2', 'yH2', 'yO2', 'yN2'}; % yH2O not in list.
    for i = 1:length(gas_list)
        gas = gas_list{i};
        if isfield(in, gas)
            gas_stream.(gas) = in.(gas)*in.n/gas_stream.n;
        end
    end
    gas_stream.yH2O = ( in.n*in.yH2O - y(1) )/gas_stream.n;

% Assign temperature unknown to gas & liq streams:
gas_stream.T = y(2);
liq_stream.T = y(2);

% Assign phases:
gas_stream.phase = "gas";
liq_stream.phase = "liq";
end