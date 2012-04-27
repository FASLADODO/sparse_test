function [ filter_avg, lres ] = avg_filter(params, rsavg)
%avg_filter Summary of this function goes here
%   Detailed explanation goes here
%   lres = avg_filter(params, rsavg)
%   lres, logical result of comparision
%   rsavg, the average result matrix
%   params, parameters to filter the result, options are as following
%       params.probe_begin, params.probe_end
%       params.time_begin, params.time_end
%       params.gamma, 
%       params.dictgamma,
%       params.lambda,
%       params.dictsize
%   {'pb1=1,pb2=93,t1=12001,t2=13500,ga=0.01,dga=1,lbd=0.01,dsize=50';}

probe_begin = 1;
probe_end = 93;

traindim = 1500;
tpoints = (1:traindim:16000);
%probe_range = 1:93:94;

%time_begin = params.time_begin;
%time_end = params.time_end;
%gamma = params.gamma;
%dictgamma = params.dictgamma;
%lambda = params.lambda;
%dictsize = params.dictsize;


%cond_key = sprintf('pb1=%d,pb2=%d,t1=%d,t2=%d,ga=%g,dga=%g,lbd=%g,dsize=%d', probe_begin, probe_end, time_begin, time_end, gamma, dictgamma, lambda, dictsize);

%compare the key of avg matrix
%cellfun(@strcmp, rsavg(:,1), repmat({con_key}, size(rsavg(:,1),1), 1));
pattern = '';

if (hasField(params, 'probe_begin') && hasField(params, 'probe_end'))
    pattern = sprintf('pb1=%d,pb2=%d,', params.probe_begin, params.probe_end);
    
else
    pattern = sprintf('pb1=%s,pb2=%s,', '.*', '.*');
    
end

if (hasField(params, 'time_begin') && hasField(params, 'time_end'))
    pt2 = sprintf('t1=%d,t2=%d,', params.time_begin, params.time_end);
    pattern = strcat(pattern, pt2);
else
    pt2 = sprintf('t1=%s,t2=%s,','.*', '.*');
    pattern = strcat(pattern, pt2);
end

if (hasField(params, 'gamma'))
    pt3 = sprintf('ga=%g,', params.gamma);
    pattern = strcat(pattern, pt3);
else
    pt3 = sprintf('ga=%s,', '.*');
    pattern = strcat(pattern, pt3);
end


if (hasField(params, 'dictgamma'))
    pt4 = sprintf('dga=%g,', params.dictgamma);
    pattern = strcat(pattern, pt4);
else
    pt4 = sprintf('dga=%s,', '.*');
    pattern = strcat(pattern, pt4);
end


if (hasField(params, 'lambda'))
    pt5 = sprintf('lbd=%g,', params.lambda);
    pattern = strcat(pattern, pt5);
else
    pt5 = sprintf('lbd=%s,', '.*');
    pattern = strcat(pattern, pt5);
end


if (hasField(params, 'dictsize'))
    pt6 = sprintf('dsize=%g', params.dictsize);
    pattern = strcat(pattern, pt6);
else
    pt6 = sprintf('dsize=%s', '.*');
    pattern = strcat(pattern, pt6);
end

disp(pattern);

%pattern = sprintf('pb1=%d,pb2=%d,t1=%d,t2=%d,ga=\w*,dga=\w*,lbd=\w*,dsize=\w*', probe_begin, probe_end, time_begin, time_end);
tmplres = cellfun(@regexp, rsavg(:,1), repmat({pattern},size(rsavg(:,1),1),1));
lres = zeros(size(tmplres));
for idx = 1:size(tmplres, 1)
    if(~isempty(tmplres(idx)) && isequal(tmplres{idx},1))
        lres(idx) = 1;
    end
        
end

 filter_avg = rsavg(logical(lres),:);


end

