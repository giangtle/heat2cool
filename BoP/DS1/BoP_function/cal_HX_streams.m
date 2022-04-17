function [hot_out, cold_out, dT_ap] = cal_HX_streams(hot_in, cold_in, dT_ap)
hot_out = hot_in; hot_out = rmfield(hot_out, "T");
cold_out = cold_in; cold_out = rmfield(cold_out, "T");
% Case 1: T approach at cold side:
hot_out.T = cold_in.T + dT_ap;
dW_hot = cal_stream_enthalpy(hot_in) - cal_stream_enthalpy(hot_out);
dW_cold = -dW_hot;
H_cold_out = cal_stream_enthalpy(cold_in) - dW_cold;
cold_out_upper = cold_in; cold_out_upper.T = hot_in.T - dT_ap;
H_cold_out_upper = cal_stream_enthalpy(cold_out_upper);
H_cold_out_lower = cal_stream_enthalpy(cold_in);
if H_cold_out > min(H_cold_out_upper,H_cold_out_lower) && H_cold_out < max(H_cold_out_upper,H_cold_out_lower)
    if abs(H_cold_out-H_cold_out_upper) >abs(H_cold_out-H_cold_out_lower)
        cold_out.T = cold_in.T; % Initial guess -> T closer to inlet
    else
        cold_out.T = cold_out_upper.T;
    end
    cold_out = cal_T_from_enthalpy(cold_out, H_cold_out);
else
    % Case 2: T approach at hot side:
    cold_out.T = hot_in.T - dT_ap;
    dW_cold = cal_stream_enthalpy(cold_in) - cal_stream_enthalpy(cold_out);
    dW_hot = -dW_cold;
    H_hot_out = cal_stream_enthalpy(hot_in) - dW_hot;
    hot_out_lower = hot_in; hot_out_lower.T = cold_in.T + dT_ap;
    H_hot_out_upper = cal_stream_enthalpy(hot_in);
    H_hot_out_lower = cal_stream_enthalpy(hot_out_lower);
    if H_hot_out > min(H_hot_out_lower, H_hot_out_upper) && H_hot_out < max(H_hot_out_lower, H_hot_out_upper)
        if abs(H_hot_out-H_hot_out_upper) > abs(H_hot_out-H_hot_out_lower)
            hot_out.T = hot_out_lower.T;
        else
            hot_out.T = hot_in.T;
        end
        hot_out = cal_T_from_enthalpy(hot_out, H_hot_out);
    else
        fprintf("dT_ap: "+dT_ap+"\n");
        fprintf("dT_ap will not be met\n");
        new_dT_ap = dT_ap + 5;
        [hot_out, cold_out] = cal_HX_streams(hot_in, cold_in, new_dT_ap);
    end
end
end