function [in, H2O_in] = cal_cooling_tower_reverse(out, T_H2O, y_H2O_in)
% 2 scenarios:  gas comes out at T > T_H2O but 100% saturated.
%               gas comes out at T = T_H2O but less than 100% saturated.
[~, liq_stream] = cal_isenthalpic_condensation(out);
if isstruct(liq_stream)
    fprintf("ERROR: input stream (in) to cooling tower is not all gaseous\n")
    return
end
% Scenario 1:
% Unknowns: 1: in.T
options = optimset("MaxFunEvals", 10000, "MaxIter", 10000, 'Display', 'off');
% Initial guess:
y0 = out.T;
% Solve:
f = @(y) equation_set(y, out, T_H2O, y_H2O_in);
[y, fval] = fsolve(f,y0,options);
if max(abs(fval))>1e-6
    return
end
[in, H2O_in] = streams_from_unknowns(y, out, T_H2O, y_H2O_in);
end

function F = equation_set(y, out, T_H2O, y_H2O_in)
[in, H2O_in] = streams_from_unknowns(y, out, T_H2O, y_H2O_in);
% 1: enthalpy(out) - enthalpy(H2O_in) - enthalpy(in)
F = cal_stream_enthalpy(in) + cal_stream_enthalpy(H2O_in) - cal_stream_enthalpy(out);
end

function [in, H2O_in] = streams_from_unknowns(y, out, T_H2O, y_H2O_in)
% H2O_in stream:
    H2O_in = struct();
    H2O_in.phase = "liq"; H2O_in.T = T_H2O; H2O_in.yH2Ol = 1;
% exit air stream:
    in = struct();
    in.phase = "gas";
    in.yH2O = y_H2O_in;
    in.n = ( 1-out.yH2O )*out.n/( 1-in.yH2O );
    % air stream compositions:
    gas_list = {'yCH4', 'yCO', 'yCO2', 'yH2', 'yO2', 'yN2'}; % yH2O not in list.
    for i = 1:length(gas_list)
        gas = gas_list{i};
        if isfield(out, gas)
            in.(gas) = out.(gas)*( 1-in.yH2O )/( 1-out.yH2O );
        end
    end
    in.T = y;
    % Also determine H2O_in.n:
    H2O_in.n = ( out.yH2O*out.n ) - ( in.yH2O*in.n );
end