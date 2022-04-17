function dP = cal_optimal_dP_Si_column(mSi, air, Dp, omega)
L=1;
dP = cal_dP_Si_column(L, mSi, air, Dp, omega);
end

function dP = cal_dP_Si_column(L, mSi, air, Dp, omega)
V = mSi/700;
Ac = V/L;
air_V = air.n*(air.T*8.31446/air.P);

mu = air_viscosity(air.T);

B = 150*(1-omega)*mu/Dp;
G = 1.75/Ac;

dP = L*(1-omega)*air_V/(Dp*omega^3*Ac)*(B+G);
end