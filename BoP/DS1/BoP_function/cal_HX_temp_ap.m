function dT_ap = cal_HX_temp_ap(hot_in, hot_out, cold_in, cold_out)
dT_ap = min(hot_in.T - cold_out.T, hot_out.T - cold_in.T);
end