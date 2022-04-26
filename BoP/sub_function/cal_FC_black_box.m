function [ex_FC, FC_black_box] = cal_FC_black_box(CH4, H2O, air, W_e, W_loss)
ex_FC = struct();
ex_FC.phase = "gas";
FC_black_box = struct();

H2O =  liq_or_vap(H2O);

ex_FC.n = CH4.n + H2O.n + air.n;
gas = "yCO2";
if isfield(air, gas)
    ex_FC.(gas) = ( air.(gas)*air.n + CH4.n )/ex_FC.n;
else
    ex_FC.(gas) = CH4.n/ex_FC.n;
end
gas = "yH2O";
if isfield(air, gas)
    ex_FC.(gas) = ( air.(gas)*air.n + H2O.n + 2*CH4.n )/ex_FC.n;
else
    ex_FC.(gas) = ( H2O.n + 2*CH4.n )/ex_FC.n;
end

ex_FC.yO2 = ( air.n*air.yO2 - 2*CH4.n )/ex_FC.n;
ex_FC.yN2 = ( air.n*air.yN2 )/ex_FC.n;

enthalpy_in = cal_stream_enthalpy(CH4) + cal_stream_enthalpy(H2O) + cal_stream_enthalpy(air);

T0 = air.T;

options = optimset("MaxFunEvals", 10000, "MaxIter", 10000, 'Display', 'off');

f = @(T) equation_set(T, ex_FC, enthalpy_in, W_e, W_loss);
[T, fval] = fsolve(f,T0,options);

if max(abs(fval))>1e-6
    return
end

ex_FC.T = T;
FC_black_box.T = T;
FC_black_box.CH4_in = CH4;
FC_black_box.H2O_in = H2O;
FC_black_box.Air_in = air;
FC_black_box.exhaust = ex_FC;
end

function F = equation_set(T, ex_FC, enthalpy_in, W_e, W_loss)
% Unknowns:     1: T
ex_FC.T = T;
F(1) = enthalpy_in - cal_stream_enthalpy(ex_FC) - W_e - W_loss;
end