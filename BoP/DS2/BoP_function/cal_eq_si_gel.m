function [y_eq, H] = cal_eq_si_gel(T, q)
kPa = 282.9425015 -1.81331662*T + -1.71093816e+02*q + 2.90522806e-03*T.^2 + 5.75477791e-01*T.*q + 8.58026192*q.^2;
y_eq = kPa/101.325;
H = y_eq/cal_water_y_sat(T);
end