function [out, reformer] = cal_reformer(in)
reformer = struct();

reformer.in = in;
% Equilibrium:
% CH4 + H2O -> CO + 3H2
% CO + H2O -> CO2 + H2

% Unknowns:     1: yCH4     2: yCO      3: yCO2     4: yH2      5:yH2O      6: T_reformer     7: n
% initial guess:
y0 = initial_guess(in);

options = optimset("MaxFunEvals", 10000, "MaxIter", 10000, 'Display', 'off');

f = @(y) full_equation_set(y, in);
[y, fval] = fsolve(f,y0,options);

if max(abs(fval))>1e-6
    y0 = second_attempt(in);
    
    options = optimset("MaxFunEvals", 5000, "MaxIter", 5000, 'Display', 'off');

    f = @(y) full_equation_set(y, in);
    [y, fval] = fsolve(f,y0,options);
    
    if max(abs(fval))>1e-6
        fprintf("Solver did not converge\n")
        return
    end
end

out = struct(); out.yCH4 = y(1); out.yCO = y(2); out.yCO2 = y(3); out.yH2 = y(4); out.yH2O = y(5); out.T = y(6); out.n = y(7);

reformer.out = out;
reformer.T = out.T;
end

function y0 = initial_guess(in)
y0 = zeros(1,7);
if in.yCH4>0
    y0(1) = in.yCH4;
else
    y0(1) = 0.01;
end
if in.yCO>0
    y0(2) = in.yCO;
else
    y0(2) = 0.01;
end
if in.yCO2>0
    y0(3) = in.yCO2;
else
    y0(3) = 0.01;
end
if in.yH2>0
    y0(4) = in.yH2;
else
    y0(4) = 0.01;
end
if in.yH2O>0
    y0(5) = in.yH2O;
else
    y0(5) = 0.01;
end
y0(6) = in.T;
y0(7) = in.n;
end

function F = full_equation_set(y, in)
% Unknowns:     1: yCH4     2: yCO      3: yCO2     4: yH2      5:yH2O      6: T_reformer     7: n
out = struct(); out.yCH4 = y(1); out.yCO = y(2); out.yCO2 = y(3); out.yH2 = y(4); out.yH2O = y(5); out.T = y(6); out.n = y(7);

K_SMR = exp(gibbs("SMR", y(6))/(-8.314472*y(6)));
K_WGS = exp(gibbs("WGS", y(6))/(-8.314472*y(6)));
% 1. CH4 + H2O -> CO + 3H2
F(1) = y(4)^3*y(2)/(y(1)*y(5)) - K_SMR;
% 2. CO + H2O -> CO2 + H2
F(2) = y(4)*y(3)/(y(2)*y(5)) - K_WGS;
% H conservation:
F(3) = y(1)*4 + y(4)*2 + y(5)*2 - (in.yCH4*4 + in.yH2*2 + in.yH2O*2)*in.n/y(7);
% C conservation:
F(4) = y(1) + y(2) + y(3) - (in.yCH4 + in.yCO + in.yCO2)*in.n/y(7);
% O conservation:
F(5) = y(2) + y(3)*2 + y(5) - (in.yCO + in.yCO2*2 + in.yH2O)*in.n/y(7);
% mole change conservation:
F(6) = ( 1-( y(1)+y(2)+y(3)+y(4)+y(5) ) ) - ( 1-( in.yCH4+in.yCO+in.yCO2+in.yH2+in.yH2O ) )*in.n/y(7);
% energy conservation:
enthalpy_in_per_mole = cal_enthalpy_per_mole(in);
enthalpy_out_per_mole = cal_enthalpy_per_mole(out);
F(7) = enthalpy_out_per_mole - enthalpy_in_per_mole*in.n/y(7);
end

function y = second_attempt(in)
y0 = initial_guess(in);

options = optimset("MaxFunEvals", 5000, "MaxIter", 5000, 'Display', 'off');

f = @(y) full_equation_set(y, in);
[y, fval] = fsolve(f,y0,options);

while max(abs(fval))>1e-6
    in.T = in.T + 50;
    
    y0 = initial_guess(in);
    
    options = optimset("MaxFunEvals", 5000, "MaxIter", 5000, 'Display', 'off');

    f = @(y) full_equation_set(y, in);
    [y, fval] = fsolve(f,y0,options);
end
end

function dH = cal_enthalpy_per_mole(in)
yCH4 = in.yCH4; yCO = in.yCO; yCO2 = in.yCO2; yH2 = in.yH2; yH2O = in.yH2O; T = in.T;
dH = ( yCH4*enthalpy("CH4",T)+yCO*enthalpy("CO",T)+yCO2*enthalpy("CO2",T)+yH2*enthalpy("H2",T)+yH2O*enthalpy("H2O",T) );
end