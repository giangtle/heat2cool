function [exhaust, FC_black_box] = cal_FC_black_box(CH4, H2O, air, W_e, W_loss)
exhaust = struct();
exhaust.phase = "gas";
FC_black_box = struct();

H2O =  water_liq_or_vap(H2O);

exhaust.n = CH4.n + H2O.n + air.n;
gas = "yCO2";
if isfield(air, gas)
    exhaust.(gas) = ( air.(gas)*air.n + CH4.n )/exhaust.n;
else
    exhaust.(gas) = CH4.n/exhaust.n;
end
gas = "yH2O";
if isfield(air, gas)
    exhaust.(gas) = ( air.(gas)*air.n + H2O.n + 2*CH4.n )/exhaust.n;
else
    exhaust.(gas) = ( H2O.n + 2*CH4.n )/exhaust.n;
end

exhaust.yO2 = ( air.n*air.yO2 - 2*CH4.n )/exhaust.n;
exhaust.yN2 = ( air.n*air.yN2 )/exhaust.n;

enthalpy_in = cal_stream_enthalpy(CH4) + cal_stream_enthalpy(H2O) + cal_stream_enthalpy(air);

T0 = air.T;

options = optimset("MaxFunEvals", 10000, "MaxIter", 10000, 'Display', 'off');

f = @(T) equation_set(T, exhaust, enthalpy_in, W_e, W_loss);
[T, fval] = fsolve(f,T0,options);

if max(abs(fval))>1e-6
    return
end

exhaust.T = T;
FC_black_box.T = T;
FC_black_box.CH4_in = CH4;
FC_black_box.H2O_in = H2O;
FC_black_box.Air_in = air;
FC_black_box.exhaust = exhaust;
end

function F = equation_set(T, exhaust, enthalpy_in, W_e, W_loss)
% Unknowns:     1: T
exhaust.T = T;
F(1) = enthalpy_in - cal_stream_enthalpy(exhaust) - W_e - W_loss;
end