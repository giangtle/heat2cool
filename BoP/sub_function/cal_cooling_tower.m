function [air_3, H2O_consumption] = cal_cooling_tower(air_4, T_env)
% Unknowns: 1: in.yH2O, H2O_in.n
options = optimset("MaxFunEvals", 10000, "MaxIter", 10000, 'Display', 'off');
% Initial guess:
y0 = [0; air_4.n*air_4.yH2O];
[air_3, H2O_consumption] = initialize_streams(T_env);
% Solve:
f = @(y) equation_set(y, air_3, H2O_consumption, air_4);
[y, fval] = fsolve(f,y0,options);
if max(abs(fval))>1e-6 || any(y<0)
    fprintf("ERROR: cal_cooling_tower.m failed to converge or return negative mole/mole fraction\n")
    air_3 = "void"; H2O_consumption = "void";
    return
end
[air_3, H2O_consumption] = streams_from_unknowns(y, air_3, H2O_consumption, air_4);
[~, liq_stream] = cal_isenthalpic_condensation(air_3);
if isstruct(liq_stream)
    fprintf("ERROR: input stream (air-3) to cooling tower is not all gaseous\n")
    air_3 = "void"; H2O_consumption = "void";
    return
end
end

function F = equation_set(y, air_3, H2O_consumption, out)
[air_3, H2O_consumption] = streams_from_unknowns(y, air_3, H2O_consumption, out);
% 1: H2O balance:
F(1) = air_3.n*air_3.yH2O + H2O_consumption.n - out.n*out.yH2O;
% 2: enthalpy(out) - enthalpy(H2O_in) - enthalpy(in)
F(2) = cal_stream_enthalpy(air_3) + cal_stream_enthalpy(H2O_consumption) - cal_stream_enthalpy(out);
end

function [air_3, H2O_consumption] = initialize_streams(T_env)
% H2O_in stream:
    H2O_consumption = struct();
    H2O_consumption.phase = "liq"; H2O_consumption.T = T_env; H2O_consumption.yH2Ol = 1;
% in stream:
    air_3 = struct();
    air_3.phase = "gas"; air_3.T = T_env+5;
end

function [air_3, H2O_consumption] = streams_from_unknowns(y, air_3, H2O_consumption, out)
% H2O_in stream:
    H2O_consumption.n = y(1);
% exit air stream:
    air_3.yH2O = y(2);
    air_3.n = ( 1-out.yH2O )*out.n/( 1-air_3.yH2O );
    % air stream compositions:
    gas_list = {'yCH4', 'yCO', 'yCO2', 'yH2', 'yO2', 'yN2'}; % yH2O not in list.
    for i = 1:length(gas_list)
        gas = gas_list{i};
        if isfield(out, gas)
            air_3.(gas) = out.(gas)*( 1-air_3.yH2O )/( 1-out.yH2O );
        end
    end
    % Also determine H2O_in.n:
    H2O_consumption.n = ( out.yH2O*out.n ) - ( air_3.yH2O*air_3.n );
end