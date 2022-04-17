function [out, mixer] = cal_mixer(in_total)
mixer = struct();

in_streams_name = fieldnames(in_total);
for k=1:numel(in_streams_name)
    mixer.(in_streams_name{k}) = in_total.(in_streams_name{k});
end


end

function out = out_stream_calculation(in_total)
in_streams_name = fieldnames(in_total);
out = struct();
sum_CH4=0; sum_CO=0; sum_CO2=0; sum_H2=0; sum_H2O=0; sum_n=0;
for k=1:numel(in_streams_name)
    sum_n   = sum_n + in_total.(in_streams_name{k}).n;
    sum_CH4 = sum_CH4 + in_total.(in_streams_name{k}).n*in_total.(in_streams_name{k}).yCH4;
    sum_CO  = sum_CO + in_total.(in_streams_name{k}).n*in_total.(in_streams_name{k}).yCO;
    
end
end