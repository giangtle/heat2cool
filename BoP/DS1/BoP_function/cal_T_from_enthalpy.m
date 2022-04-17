function in = cal_T_from_enthalpy(in, H)
% Initial guess:
if isfield(in, "T")
    T0 = in.T;
else
    T0 = 200+273.15;
end
% Solve:
options = optimset("MaxFunEvals", 10000, "MaxIter", 10000, 'Display', 'off');
f = @(T) eqn(T, H, in);
[T, fval] = fsolve(f, T0, options);
if max(abs(fval))>1e-4
    return
end
in.T = T;
    function F = eqn(T, H, in)
        in.T = T;
        F = cal_stream_enthalpy(in) - H;
    end
end