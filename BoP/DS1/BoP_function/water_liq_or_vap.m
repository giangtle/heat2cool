function out = water_liq_or_vap(H2O_in)
out = struct();
% check if "phase" field is in H2O_in:
if isfield(H2O_in, "phase")
    if H2O_in.phase == "liq"
        if H2O_in.yH2Ol ~= 1
            fprintf("Not pure water");
            return
        end
    elseif H2O_in.phase == "gas"
        if H2O_in.yH2O ~= 1
            fprintf("Not pure water");
            return
        end
    else
        fprintf("Wrong identifier for phase")
    end
end
% Find correct phase label:
out.n = H2O_in.n;
out.T = H2O_in.T;
if H2O_in.T > (100+273.15)
    out.yH2O = 1;
    out.phase = "gas";
else
    out.yH2Ol = 1;
    out.phase = "liq";
end
end