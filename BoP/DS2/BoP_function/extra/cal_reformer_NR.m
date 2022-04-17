function [out, reformer] = cal_reformer_NR(in)
reformer = struct();

reformer.in = in;
% Equilibrium:
% CH4 + H2O -> CO + 3H2
% CO + H2O -> CO2 + H2

% Unknowns:     1: yCH4     2: yCO      3: yCO2     4: yH2      5:yH2O      6: T_reformer     7: n
dH = zeros(size(in.T:-25:in.T-400));
pos = 1;
y = [0.0005; 0.2; 0.2; 0.2; 0.2; in.T; in.n];

% for T=in.T:-25:in.T-400
%     y(6) = T;
%     [F, Jacobian]= find_NR_Jacobian(y, reformer);
%     error=max(abs(F));
%     iteration = 1;
%     while (error>1e-7 && iteration<15)
%         y = y - Jacobian\F;
%         [F, Jacobian]= find_NR_Jacobian(y, reformer);
%         error=max(abs(F));
%         iteration = iteration+1;
%     end
%     out = struct(); out.yCH4 = y(1); out.yCO = y(2); out.yCO2 = y(3); out.yH2 = y(4); out.yH2O = y(5); out.T = y(6); out.n = y(7);
%     enthalpy_in_per_mole = cal_enthalpy(in);
%     enthalpy_out_per_mole = cal_enthalpy(out);
%     dH(pos) = enthalpy_out_per_mole - enthalpy_in_per_mole*in.n/out.n;
%     pos = pos+1;
% end

    [F, Jacobian]= find_NR_Jacobian(y, reformer);
    error=max(abs(F));
    iteration = 1;
    while (error>1e-7 && iteration<15)
        y = y - Jacobian\F;
        [F, Jacobian]= find_NR_Jacobian(y, reformer);
        error=max(abs(F));
        iteration = iteration+1;
    end

    [F, Jacobian]= find_NR_Jacobian_2(y, reformer);
    error=max(abs(F));
    iteration = 1;
    while (error>1e-7 && iteration<50)
        y = y - Jacobian\F;
        [F, Jacobian]= find_NR_Jacobian_2(y, reformer);
        error=max(abs(F));
        iteration = iteration+1;
    end
    
if error > 1e-6
    y0 = y;
    options = optimset("MaxFunEvals", 10000, "MaxIter", 10000, 'Display', 'off');
    f = @(y) equation_set(y, reformer);
    [y, fval] = fsolve(f,y0,options);
end

out = struct(); out.yCH4 = y(1); out.yCO = y(2); out.yCO2 = y(3); out.yH2 = y(4); out.yH2O = y(5); out.T = y(6); out.n = y(7);

reformer.out = out;
reformer.T = out.T;
end

function [F, Jacobian]= find_NR_Jacobian(y, reformer)
size=length(y); F=zeros(size,1); Jacobian=eye(size,size);

in = reformer.in;
out = struct(); out.yCH4 = y(1); out.yCO = y(2); out.yCO2 = y(3); out.yH2 = y(4); out.yH2O = y(5); out.T = y(6); out.n = y(7);
% Equilibrium:
% 1. CH4 + H2O -> CO + 3H2
% 2. CO + H2O -> CO2 + H2

dT = 0.001;
K_SMR = exp(gibbs("SMR", out.T)/(-8.314472*out.T));
K_WGS = exp(gibbs("WGS", out.T)/(-8.314472*out.T));

for eqncount=1
    % 1. CH4 + H2O -> CO + 3H2
    F(eqncount) = out.yH2^3*out.yCO/(out.yCH4*out.yH2O) - K_SMR;
    Jacobianrow = zeros(1, size);
    Jacobianrow(1) = -out.yH2^3*out.yCO/(out.yCH4^2*out.yH2O);
    Jacobianrow(2) = out.yH2^3/(out.yCH4*out.yH2O);
    Jacobianrow(4) = 3*out.yH2^2*out.yCO/(out.yCH4*out.yH2O);
    Jacobianrow(5) = -out.yH2^3*out.yCO/(out.yCH4*out.yH2O^2);
    Jacobian(eqncount, :)   = Jacobianrow;
end

for eqncount=2
    % 2. CO + H2O -> CO2 + H2
    F(eqncount) = out.yH2*out.yCO2/(out.yCO*out.yH2O) - K_WGS;
    Jacobianrow = zeros(1, size);
    Jacobianrow(2) = -out.yH2*out.yCO2/(out.yCO^2*out.yH2O);
    Jacobianrow(3) = out.yH2/(out.yCO*out.yH2O);
    Jacobianrow(4) = out.yCO2/(out.yCO*out.yH2O);
    Jacobianrow(5) = -out.yH2*out.yCO2/(out.yCO*out.yH2O^2);
    Jacobian(eqncount, :)   = Jacobianrow;
end

for eqncount = 3
    % H conservation:
    F(eqncount) = out.yCH4*4 + out.yH2*2 + out.yH2O*2 - (in.yCH4*4 + in.yH2*2 + in.yH2O*2)*in.n/out.n;
    Jacobianrow = zeros(1, size);
    Jacobianrow(1) = 4;
    Jacobianrow(4) = 2;
    Jacobianrow(5) = 2;
    Jacobianrow(7) = (in.yCH4*4 + in.yH2*2 + in.yH2O*2)*in.n/out.n^2;
    Jacobian(eqncount, :)   = Jacobianrow;
end

for eqncount = 4
    % C conservation:
    F(eqncount) = out.yCH4 + out.yCO + out.yCO2 - (in.yCH4 + in.yCO + in.yCO2)*in.n/out.n;
    Jacobianrow = zeros(1, size);
    Jacobianrow(1) = 1;
    Jacobianrow(2) = 1;
    Jacobianrow(3) = 1;
    Jacobianrow(7) = (in.yCH4 + in.yCO + in.yCO2)*in.n/out.n^2;
    Jacobian(eqncount, :)   = Jacobianrow;
end

for eqncount = 5
    % O conservation:
    F(eqncount) = out.yCO + out.yCO2*2 + out.yH2O - (in.yCO + in.yCO2*2 + in.yH2O)*in.n/out.n;
    Jacobianrow = zeros(1, size);
    Jacobianrow(2) = 1;
    Jacobianrow(3) = 2;
    Jacobianrow(5) = 1;
    Jacobianrow(7) = (in.yCO + in.yCO2*2 + in.yH2O)*in.n/out.n^2;
    Jacobian(eqncount, :)   = Jacobianrow;
end

% for eqncount = 6
%     % adiabatic condition:
%     enthalpy_in_per_mole = cal_enthalpy(in);
%     enthalpy_out_per_mole = cal_enthalpy(out);
%     F(eqncount) = enthalpy_out_per_mole - enthalpy_in_per_mole*in.n/out.n;
%     Jacobianrow = zeros(1, size);
%     Jacobianrow(1) = enthalpy("CH4",out.T);
%     Jacobianrow(2) = enthalpy("CO",out.T);
%     Jacobianrow(3) = enthalpy("CO2",out.T);
%     Jacobianrow(4) = enthalpy("H2",out.T);
%     Jacobianrow(5) = enthalpy("H2O",out.T);
%     dT_out = out; dT_out.T = out.T+dT;
%     Jacobianrow(6) = (cal_enthalpy(dT_out)-cal_enthalpy(out))/dT;
%     Jacobianrow(7) = enthalpy_in_per_mole*in.n/out.n^2;
%     Jacobian(eqncount, :)   = Jacobianrow;
% end

for eqncount = 7
    % mole fraction added to 1:
    F(eqncount) = ( 1-( out.yCH4+out.yCO+out.yCO2+out.yH2+out.yH2O ) ) - ( 1-( in.yCH4+in.yCO+in.yCO2+in.yH2+in.yH2O ) )*in.n/out.n;
    Jacobianrow = zeros(1, size);
    Jacobianrow(1) = -1;
    Jacobianrow(2) = -1;
    Jacobianrow(3) = -1;
    Jacobianrow(4) = -1;
    Jacobianrow(5) = -1;
    Jacobianrow(7) = ( 1-( in.yCH4+in.yCO+in.yCO2+in.yH2+in.yH2O ) )*in.n/out.n^2;
    Jacobian(eqncount, :)   = Jacobianrow;
end
end

function [F, Jacobian]= find_NR_Jacobian_2(y, reformer)
size=length(y); F=zeros(size,1); Jacobian=eye(size,size);

in = reformer.in;
out = struct(); out.yCH4 = y(1); out.yCO = y(2); out.yCO2 = y(3); out.yH2 = y(4); out.yH2O = y(5); out.T = y(6); out.n = y(7);
% Equilibrium:
% 1. CH4 + H2O -> CO + 3H2
% 2. CO + H2O -> CO2 + H2

dT = 0.001;
K_SMR = exp(gibbs("SMR", out.T)/(-8.314472*out.T));
K_WGS = exp(gibbs("WGS", out.T)/(-8.314472*out.T));

for eqncount=1
    % 1. CH4 + H2O -> CO + 3H2
    F(eqncount) = out.yH2^3*out.yCO/(out.yCH4*out.yH2O) - K_SMR;
    Jacobianrow = zeros(1, size);
    Jacobianrow(1) = -out.yH2^3*out.yCO/(out.yCH4^2*out.yH2O);
    Jacobianrow(2) = out.yH2^3/(out.yCH4*out.yH2O);
    Jacobianrow(4) = 3*out.yH2^2*out.yCO/(out.yCH4*out.yH2O);
    Jacobianrow(5) = -out.yH2^3*out.yCO/(out.yCH4*out.yH2O^2);
    Jacobianrow(6) = -( exp(gibbs("SMR", (out.T+dT))/(-8.314472*(out.T+dT)))-exp(gibbs("SMR", out.T)/(-8.314472*out.T)) )/dT;
    Jacobian(eqncount, :)   = Jacobianrow;
end

for eqncount=2
    % 2. CO + H2O -> CO2 + H2
    F(eqncount) = out.yH2*out.yCO2/(out.yCO*out.yH2O) - K_WGS;
    Jacobianrow = zeros(1, size);
    Jacobianrow(2) = -out.yH2*out.yCO2/(out.yCO^2*out.yH2O);
    Jacobianrow(3) = out.yH2/(out.yCO*out.yH2O);
    Jacobianrow(4) = out.yCO2/(out.yCO*out.yH2O);
    Jacobianrow(5) = -out.yH2*out.yCO2/(out.yCO*out.yH2O^2);
    Jacobianrow(6) = -( exp(gibbs("WGS", (out.T+dT))/(-8.314472*(out.T+dT)))-exp(gibbs("WGS", out.T)/(-8.314472*out.T)) )/dT;
    Jacobian(eqncount, :)   = Jacobianrow;
end

for eqncount = 3
    % H conservation:
    F(eqncount) = out.yCH4*4 + out.yH2*2 + out.yH2O*2 - (in.yCH4*4 + in.yH2*2 + in.yH2O*2)*in.n/out.n;
    Jacobianrow = zeros(1, size);
    Jacobianrow(1) = 4;
    Jacobianrow(4) = 2;
    Jacobianrow(5) = 2;
    Jacobianrow(7) = (in.yCH4*4 + in.yH2*2 + in.yH2O*2)*in.n/out.n^2;
    Jacobian(eqncount, :)   = Jacobianrow;
end

for eqncount = 4
    % C conservation:
    F(eqncount) = out.yCH4 + out.yCO + out.yCO2 - (in.yCH4 + in.yCO + in.yCO2)*in.n/out.n;
    Jacobianrow = zeros(1, size);
    Jacobianrow(1) = 1;
    Jacobianrow(2) = 1;
    Jacobianrow(3) = 1;
    Jacobianrow(7) = (in.yCH4 + in.yCO + in.yCO2)*in.n/out.n^2;
    Jacobian(eqncount, :)   = Jacobianrow;
end

for eqncount = 5
    % O conservation:
    F(eqncount) = out.yCO + out.yCO2*2 + out.yH2O - (in.yCO + in.yCO2*2 + in.yH2O)*in.n/out.n;
    Jacobianrow = zeros(1, size);
    Jacobianrow(2) = 1;
    Jacobianrow(3) = 2;
    Jacobianrow(5) = 1;
    Jacobianrow(7) = (in.yCO + in.yCO2*2 + in.yH2O)*in.n/out.n^2;
    Jacobian(eqncount, :)   = Jacobianrow;
end

for eqncount = 6
    % adiabatic condition:
    enthalpy_in_per_mole = cal_enthalpy(in);
    enthalpy_out_per_mole = cal_enthalpy(out);
    F(eqncount) = enthalpy_out_per_mole - enthalpy_in_per_mole*in.n/out.n;
    Jacobianrow = zeros(1, size);
    dT_out = out; dT_out.T = out.T+dT;
    Jacobianrow(6) = (cal_enthalpy(dT_out)-cal_enthalpy(out))/dT*20;
    Jacobian(eqncount, :)   = Jacobianrow;
end

for eqncount = 7
    % mole fraction added to 1:
    F(eqncount) = ( 1-( out.yCH4+out.yCO+out.yCO2+out.yH2+out.yH2O ) ) - ( 1-( in.yCH4+in.yCO+in.yCO2+in.yH2+in.yH2O ) )*in.n/out.n;
    Jacobianrow = zeros(1, size);
    Jacobianrow(1) = -1;
    Jacobianrow(2) = -1;
    Jacobianrow(3) = -1;
    Jacobianrow(4) = -1;
    Jacobianrow(5) = -1;
    Jacobianrow(7) = ( 1-( in.yCH4+in.yCO+in.yCO2+in.yH2+in.yH2O ) )*in.n/out.n^2;
    Jacobian(eqncount, :)   = Jacobianrow;
end
end

function F = equation_set(y, reformer)
% Unknowns:     1: yCH4     2: yCO      3: yCO2     4: yH2      5:yH2O      6: T_reformer     7: n
in = reformer.in;
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
enthalpy_in_per_mole = cal_enthalpy(in);
enthalpy_out_per_mole = cal_enthalpy(out);
F(7) = enthalpy_out_per_mole - enthalpy_in_per_mole*in.n/y(7);
end

function dH = cal_enthalpy(in)
yCH4 = in.yCH4; yCO = in.yCO; yCO2 = in.yCO2; yH2 = in.yH2; yH2O = in.yH2O; T = in.T;
dH = ( yCH4*enthalpy("CH4",T)+yCO*enthalpy("CO",T)+yCO2*enthalpy("CO2",T)+yH2*enthalpy("H2",T)+yH2O*enthalpy("H2O",T) );
end