function dH = cal_stream_enthalpy(in)
dH = 0;
gas_list = {'yCH4', 'yCO', 'yCO2', 'yH2', 'yH2O', 'yO2', 'yN2'};
liq_list = {'yH2Ol'};
if in.phase == "gas"
    for i = 1:length(gas_list)
        gas = gas_list{i};
        if isfield(in, gas)
            dH = dH + in.(gas)*in.n*enthalpy(gas(2:end),in.T);
        end
    end
end
if in.phase == "liq"
    for i = 1:length(liq_list)
        liq = liq_list{i};
        if isfield(in, liq)
            dH = dH + in.(liq)*in.n*enthalpy(liq(2:end),in.T);
        end
    end
end
if in.phase == "liq-gas"
    for i = 1:length(gas_list)
        gas = gas_list{i};
        if isfield(in, gas)
            dH = dH + in.(gas)*in.n_gas*enthalpy(gas(2:end),in.T);
        end
    end
    for i = 1:length(liq_list)
        liq = liq_list{i};
        if isfield(in, liq)
            dH = dH + in.(liq)*in.n_liq*enthalpy(liq(2:end),in.T);
        end
    end
end
end