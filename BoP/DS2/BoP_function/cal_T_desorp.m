function [T_desorp, ex_out] = cal_T_desorp(ex_in, T_desorp)
y_EQ = solve_y_EQ(ex_in, T_desorp);
% If y_EQ is greater than the saturated water mole fraction, increase desorption temperature.
while y_EQ > cal_water_y_sat(T_desorp)
    T_desorp = T_desorp+10;
    y_EQ = solve_y_EQ(ex_in, T_desorp);
end
[~, ex_out]=equation_set(y_EQ, ex_in, T_desorp);
end

function y_EQ = solve_y_EQ(ex_in, T_desorp)
% Unknowns: 1: ex_out.yH2O/ y_EQ
options = optimset("MaxFunEvals", 10000, "MaxIter", 10000, 'Display', 'off');
% Initial guess:
y0 = ex_in.yH2O*1.1;
% Solve:
f = @(y) equation_set(y, ex_in, T_desorp);
[y_EQ, fval] = fsolve(f,y0,options);
if max(abs(fval))>1e-6
    return
end
end

function ex_out = cal_exhaust_out(ex_in, y_EQ, T_desorp)
% ex_out stream:
ex_out = ex_in;
ex_out.T = T_desorp;
ex_out.yH2O = y_EQ;
ex_out.n = ex_in.n*( 1-ex_in.yH2O )/( 1-ex_out.yH2O );
% Other component mole fraction:
gas_list = {'yCH4', 'yCO', 'yCO2', 'yH2', 'yO2', 'yN2'}; % yH2O not in list.
for i = 1:length(gas_list)
    gas = gas_list{i};
    if isfield(ex_in, gas)
        ex_out.(gas) = ex_in.(gas)*ex_in.n/ex_out.n;
    end
end  
end

function [F, ex_out]=equation_set(y, ex_in, T_desorp)
% water MW:
H2O_MW = 18.01528;  % [g/mol]
% Heat of adsorption:
Q_adsorption = 2800;    % [J/g]
% delta enthalpy:
ex_out = cal_exhaust_out(ex_in, y, T_desorp);
dH = cal_stream_enthalpy(only_air(ex_in))-cal_stream_enthalpy(only_air(ex_out));
% delta water mass:
dmH2O_desorption = (ex_in.n*ex_in.yH2O-ex_out.n*ex_out.yH2O)*H2O_MW*Q_adsorption;
% Energy balance:
F = dH + dmH2O_desorption;
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