function [gas_stream, liq_stream] = cal_isothermal_condensation(in)
if in.phase == "liq"
    gas_stream = 0;
    liq_stream = in;
else
    yH2O_sat = cal_water_y_sat(in.T);
    if in.yH2O <= yH2O_sat
        gas_stream = in;
        liq_stream = 0;
    else
        [gas_stream, liq_stream] = isothermal_condesation(in);
    end
end
end

function [gas_stream, liq_stream] = isothermal_condesation(in)
y_sat = cal_water_y_sat(in.T);
n_gas_stream = in.n*(1-in.yH2O)/(1-y_sat);
n_liq_stream = in.n- n_gas_stream;
% Assign streams:
% liq stream:
liq_stream = struct();
liq_stream.phase = "liq";
liq_stream.T = in.T;
liq_stream.n = n_liq_stream;
liq_stream.yH2Ol = 1;
% gas stream:
gas_stream = struct();
gas_stream.phase = "gas";
gas_stream.T = in.T;
gas_stream.n = n_gas_stream;
gas_stream.yH2O = y_sat;
gas_list = {'yCH4', 'yCO', 'yCO2', 'yH2', 'yO2', 'yN2'}; % yH2O not in list.
for i = 1:length(gas_list)
    gas = gas_list{i};
    if isfield(in, gas)
        gas_stream.(gas) = in.(gas)*in.n/gas_stream.n;
    end
end
end