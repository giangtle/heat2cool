function [out, H2O_in] = cal_cooling_tower(in, T_H2O)
% 2 scenarios:  gas comes out at T > T_H2O but 100% saturated.
%               gas comes out at T = T_H2O but less than 100% saturated.
[~, liq_stream] = cal_condensation(in);
if isstruct(liq_stream)
    fprintf("ERROR: input stream (in) to cooling tower is not all gaseous\n")
    return
end
if T_H2O > in.T
    fprintf("ERROR: water temperature higher than air, no cooling\n")
    return
end
% Scenario 1:
% Unknowns: 1: H2O_in.n     2: out.T
options = optimset("MaxFunEvals", 10000, "MaxIter", 10000, 'Display', 'off');
% Initial guess:
y0 = [0, in.T];
% Solve:
f = @(y) equation_set_1(y, in, T_H2O);
[y, fval] = fsolve(f,y0,options);
if max(abs(fval))>1e-6
    return
end
if y(2) > T_H2O
    [out, H2O_in] = streams_from_unknowns_1(y, in, T_H2O);
    return
end

% Scenario 2:
% Unknowns: 1: H2O_in.n
% Initial guess:
y0 = 0;
% Solve:
f = @(y) equation_set_2(y, in, T_H2O);
[y, fval] = fsolve(f,y0,options);
if max(abs(fval))>1e-6
    return
end
[out, H2O_in] = streams_from_unknowns_2(y, in, T_H2O);
end

function F = equation_set_1(y, in, T_H2O)
[out, H2O_in] = streams_from_unknowns_1(y, in, T_H2O);
% Equations:
if isfield(in, "n")
    % 1: yH2O = yH2O_sat = exp(-gibbs("evap", y(2))/(8.314472*y(2))
    F(1) = ( in.n*in.yH2O + y(1) )/( in.n + y(1) ) - cal_water_y_sat(y(2));
    % 2: enthalpy(out) - enthalpy(H2O_in) - enthalpy(in)
    F(2) = cal_stream_enthalpy(out) - cal_stream_enthalpy(H2O_in) - cal_stream_enthalpy(in);
else
    % 1: yH2O = yH2O_sat = exp(-gibbs("evap", y(2))/(8.314472*y(2))
    F(1) = ( in.yH2O + y(1) )/( 1 + y(1) ) - cal_water_y_sat(y(2));
    % 2: enthalpy(out) - enthalpy(H2O_in) - enthalpy(in)
    in.n = 1;
    out.n = 1+y(1);
    F(2) = cal_stream_enthalpy(out) - cal_stream_enthalpy(H2O_in) - cal_stream_enthalpy(in);
end
end

function F = equation_set_2(y, in, T_H2O)
[out, H2O_in] = streams_from_unknowns_2(y, in, T_H2O);
if isfield(in, "n")
    % 1: enthalpy(out) - enthalpy(H2O_in) - enthalpy(in)
    F = cal_stream_enthalpy(out) - cal_stream_enthalpy(H2O_in) - cal_stream_enthalpy(in);
else
    in.n = 1;
    out.n = 1+y(1);
    F = cal_stream_enthalpy(out) - cal_stream_enthalpy(H2O_in) - cal_stream_enthalpy(in);
end
end

function [out, H2O_in] = streams_from_unknowns_1(y, in, T_H2O)
% H2O_in stream:
    H2O_in = struct();
    H2O_in.phase = "liq"; H2O_in.n = y(1); H2O_in.T = T_H2O; H2O_in.yH2Ol = 1;
% exit air stream:
    out = struct();
    out.phase = "gas";
    if isfield(in, "n")
        out.n = in.n + y(1);
        % air stream compositions:
        gas_list = {'yCH4', 'yCO', 'yCO2', 'yH2', 'yO2', 'yN2'}; % yH2O not in list.
        for i = 1:length(gas_list)
            gas = gas_list{i};
            if isfield(in, gas)
                out.(gas) = in.(gas)*in.n/out.n;
            end
        end
        out.yH2O = ( in.yH2O*in.n + y(1) )/out.n;
    else
        % air stream compositions:
        gas_list = {'yCH4', 'yCO', 'yCO2', 'yH2', 'yO2', 'yN2'}; % yH2O not in list.
        for i = 1:length(gas_list)
            gas = gas_list{i};
            if isfield(in, gas)
                out.(gas) = in.(gas)*1/( 1+y(1) );
            end
        end
        out.yH2O = ( in.yH2O + y(1) )/( 1+y(1) );
    end
    out.T = y(2);
end

function [out, H2O_in] = streams_from_unknowns_2(y, in, T_H2O)
% H2O_in stream:
    H2O_in = struct();
    H2O_in.phase = "liq"; H2O_in.n = y; H2O_in.T = T_H2O; H2O_in.yH2Ol = 1;
% exit air stream:
    out = struct();
    out.phase = "gas";
    if isfield(in, "n")
        out.n = in.n + y;
        % air stream compositions:
        gas_list = {'yCH4', 'yCO', 'yCO2', 'yH2', 'yO2', 'yN2'}; % yH2O not in list.
        for i = 1:length(gas_list)
            gas = gas_list{i};
            if isfield(in, gas)
                out.(gas) = in.(gas)*in.n/out.n;
            end
        end
        out.yH2O = ( in.yH2O*in.n + y )/out.n;
        out.T = T_H2O;
    else
        % air stream compositions:
        gas_list = {'yCH4', 'yCO', 'yCO2', 'yH2', 'yO2', 'yN2'}; % yH2O not in list.
        for i = 1:length(gas_list)
            gas = gas_list{i};
            if isfield(in, gas)
                out.(gas) = in.(gas)/( 1+y );
            end
        end
        out.yH2O = ( in.yH2O + y )/( 1+y );
        out.T = T_H2O;
    end
end