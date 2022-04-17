function [dP, power, R_column, t_res] = cal_dP_Si_column(L, mSi, air, Dp, omega)
V = mSi/700;
Ac = V/L;
R_column = sqrt(Ac/pi);
air_V = air.n*(air.T*8.31446/air.P);

mu = air_viscosity(air.T);
ro = air_density(air);

vs = air_V/Ac;

Rep = ro*vs*Dp/(1-omega)/mu;

fp = 150/Rep + 1.75;

dP = fp*L*ro*vs^2*(1-omega)/(Dp*omega^3);
power = dP*air_V;

dP = dP/101325;
t_res = V/(air_V);
end