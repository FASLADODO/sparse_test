function [probe_params, resmap]=data_proc(regstr)
% usage, input the filename.mat files you want to parse, and get return as
% the following format
% format of probe_params
% <data_start_time, data_end_time, gamma, lambda, dictsize, inner_index, train_err, test_err>
% this function may need change when there are more than one probe_ranges,
% e.g. [1, 100], [101, 200], ... etc.

matfiles = dir(regstr);
tol_matrix = [];

for i = 1:numel(matfiles)
    outstr = sprintf('%s\n', matfiles(i).name);
    disp(outstr);
   
    clearvars -except matfiles i tol_matrix;
    loadstrc = load(matfiles(i).name, 'objs');
    res_cells = loadstrc.objs;
    [~, ~, ~, rmatrix, resmap]=res_extract(res_cells);
    tol_matrix = [tol_matrix; rmatrix];
    
end

rmatrix = tol_matrix;


% <itest, probe_range_begin, probe_range_end, time_range_begin, time_range_end, gamma, dictgamma, lambda, dictsize, iter_num, step1_obj, step2_obj, train_err, test_err>
% for each dataset compute the 
% find the unique probe range
%probes = unique(rmatrix(:, 2:5), 'rows');

% values from 6 to 9 are <gamma, dgamma, lambda, dictsize, inner_iter_num>
%params = unique(rmatrix(:, 6:9), 'rows');

propars = unique(rmatrix(:, 2:10), 'rows');
probe_params = zeros(size(propars, 1), 9);

for ipms = 1:size(propars, 1)
    pstart_idx = rmatrix(:,2) == propars(ipms,1);
    pend_idx = rmatrix(:,3) == propars(ipms, 2);
    tstart_idx = rmatrix(:, 4) == propars(ipms, 3);
    tend_idx = rmatrix(:, 5) == propars(ipms, 4);
    probe_idx = pstart_idx .* pend_idx .* tstart_idx .* tend_idx;
    
    gamma_idx = rmatrix(:, 6) == propars(ipms, 5);
    dgamma_idx = rmatrix(:, 7) == propars(ipms, 6);
    lambda_idx = rmatrix(:, 8) == propars(ipms, 7);
    dict_idx = rmatrix(:, 9) == propars(ipms, 8);
    inner_idx = rmatrix(:, 10) == propars(ipms, 9);
    params_idx = gamma_idx .* dgamma_idx .* lambda_idx .* dict_idx .* inner_idx;
    all_idx = probe_idx .* params_idx;
    mtrain_err = mean(rmatrix(logical(all_idx), 13));
    mtest_err = mean(rmatrix(logical(all_idx), 14));
    
    probe_params(ipms, 1) = propars(ipms, 3); % record the start time
    probe_params(ipms, 2) = propars(ipms, 4); % record the end time
    probe_params(ipms, 3) = propars(ipms, 5); % gamma
    probe_params(ipms, 4) = propars(ipms, 6); % dictgamma
    probe_params(ipms, 5) = propars(ipms, 7); % lambda
    probe_params(ipms, 6) = propars(ipms, 8); % dictsize
    probe_params(ipms, 7) = propars(ipms, 9); % inner iteration
    
    probe_params(ipms, 8) = mtrain_err; % train err
    probe_params(ipms, 9) = mtest_err; % test err
    
    
end


%{
% store the statistical result group by probe and params
probe_params = zeros(size(probes,1)*size(params,1)*0.5, 8);

for iprobe = 1:(size(probes,1))
    
    pstart_idx = rmatrix(:,2) == probes(iprobe,1);
    pend_idx = rmatrix(:,3) == probes(iprobe, 2);
    tstart_idx = rmatrix(:, 4) == probes(iprobe, 3);
    tend_idx = rmatrix(:, 5) == probes(iprobe, 4);
    probe_idx = pstart_idx .* pend_idx .* tstart_idx .* tend_idx;
    
   
    for iparam = 1:size(params,1)
        gamma_idx = rmatrix(:, 6) == params(iparam, 1);
        lambda_idx = rmatrix(:, 7) == params(iparam, 2);
        dict_idx = rmatrix(:, 8) == params(iparam, 3);
        inner_idx = rmatrix(:, 9) == params(iparam, 4);
        params_idx = gamma_idx .* lambda_idx .* dict_idx .* inner_idx;
        all_idx = probe_idx .* params_idx;
        mtrain_err = mean(rmatrix(logical(all_idx), 12));
        mtest_err = mean(rmatrix(logical(all_idx), 13));
        
        probe_params((iprobe-1)*size(probes,1)+iparam, 1) = probes(iprobe, 3); % record the start time
        probe_params((iprobe-1)*size(probes,1)+iparam, 2) = probes(iprobe, 4);
        probe_params((iprobe-1)*size(probes,1)+iparam, 3) = params(iparam, 1);
        probe_params((iprobe-1)*size(probes,1)+iparam, 4) = params(iparam, 2);
        probe_params((iprobe-1)*size(probes,1)+iparam, 5) = params(iparam, 3);
        probe_params((iprobe-1)*size(probes,1)+iparam, 6) = params(iparam, 4);
        probe_params((iprobe-1)*size(probes,1)+iparam, 7) = mtrain_err;
        probe_params((iprobe-1)*size(probes,1)+iparam, 8) = mtest_err;
    end

end
%}
        


end