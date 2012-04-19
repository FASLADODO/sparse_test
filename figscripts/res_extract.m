function [params_cells, iteration_cells, dataset_cells, res_matrix, res_map]=res_extract(res_cells)
% res_matrix, with the format <itest, probe_range_begin, probe_range_end, time_range_begin, time_range_end, gamma, lambda, dictsize, iter_num, step1_obj, step2_obj, train_err, test_err>
% dataset_cells, for a given dataset, get all its results
% iteration_cells, group results by iteration
% params_cells, group results by parameters
% res_map, the map contains the learned dictionary and leanred cofficients
% matrix, with all the parameters as the key
%
% res_cell, is a 6 dimension cell matrix storing all the temporal result
% objs{itest, kk, jj, igamma, ilambda, idict} = {params, obj_func, outD, outX}; //line 80, sparse_test.m

% $itest, is the current iteration number
% $kk, is the num of kth row range of  observed matrix Y, the original observed Y (observed probes * time_slices) is too large for computing, thus we divide it into sub range, $kk is the mark for set of observed probes
% $jj, is the column range of observed matrix Y, i.e. the time dimension
% $igamma, is the ith gamma value for 1-step optimization
% $ilambda, is the ith lambda value for 2-step optimization
% $idict, is the ith dictionary size value

% params, is a struct which keeps all the input params, including $itest, $kk, $jj, $igamma, $ilambda, $idict,$initDict, etc...
% obj_func, is a <4 \times params.iternum> matrix, each column contains the test result of a iteration. <1, :> is objective function of 1-step feature sign, 
% <2, :> is objective function of 2-step feature sign, <3, :> is rmse for train set, <4, :> is rmse for test set.
% outD, is the learned dictionary
% outX, is the learned coefficients matrix 

[test_iter, probe_range, time_range, gamma_range, dgamma_range, lambda_range, dict_range] = size(res_cells);


% for a given parameters pair, get all its test results on different
% datasets
%params_cells = cell(size(gamma_range), size(lambda_range), size(dict_range));
params_cells = cell(gamma_range, dgamma_range, lambda_range, dict_range);
for igamma = 1:gamma_range
    for idgamma =1:dgamma_range
        for ilambda = 1:lambda_range
            for idict = 1:dict_range
                params_cells{igamma, idgamma, ilambda, idict} = res_cells(:, :, :, igamma, idgamma, ilambda, idict);
            end
        end
    end
end


% group result by iteration
iteration_cells = reshape(res_cells, test_iter, numel(res_cells(1,:)));


% for a given dataset, get all the results
dataset_cells = cell(probe_range, time_range);
for iprobe = 1:probe_range
    for itime = 1:time_range
        dataset_cells{iprobe, itime} = res_cells(:, iprobe, itime, :, :, :, :);
    end
end

% extract all the result into a matrix
% [itest, probe_range_begin, probe_range_end, time_range_begin, time_range_end, gamma, dgamma, lambda, dictsize, iter_num, step1_obj, step2_obj, train_err, test_err]
res_matrix = zeros(numel(res_cells), 14);
res_map = containers.Map();

num_rec = 1;
for itest = 1:test_iter
    for iprobe = 1:probe_range
        for itime = 1:time_range
            for igamma = 1:gamma_range
                for idgamma = 1:dgamma_range
                    for ilambda = 1:lambda_range
                        for idict = 1:dict_range
                            if ~isempty(res_cells{itest, iprobe, itime, igamma, idgamma, ilambda, idict})

                                inner_iter = res_cells{itest, iprobe, itime, igamma, idgamma, ilambda, idict}{1}.iternum; % get the num of inner iterations, i.e. params.iternum

                                cgamma = res_cells{itest, iprobe, itime, igamma, idgamma, ilambda, idict}{1}.gamma;
                                cdgamma = res_cells{itest, iprobe, itime, igamma, idgamma, ilambda, idict}{1}.dictgamma;
                                clambda = res_cells{itest, iprobe, itime, igamma, idgamma, ilambda, idict}{1}.lambda;
                                cdict = size(res_cells{itest, iprobe, itime, igamma, idgamma, ilambda, idict}{1}.initDict, 2);
                                mapkey = sprintf('iter=%d,prb=%d,t=%d,ga=%g,dga=%g,la=%g,dict=%d', itest, iprobe, itime, cgamma, cdgamma, clambda, cdict);
                                outD = res_cells{itest, iprobe, itime, igamma, idgamma, ilambda, idict}{3};
                                outX = res_cells{itest, iprobe, itime, igamma, idgamma, ilambda, idict}{4};
                                res_map(mapkey) = {outD, outX};


                                for iiter = 1:inner_iter
                                    res_matrix(num_rec, 1) = itest;

                                    prange = size(res_cells{itest, iprobe, itime, igamma, idgamma, ilambda, idict}{1}.probe, 2);
                                    res_matrix(num_rec, 2) = res_cells{itest, iprobe, itime, igamma, idgamma, ilambda, idict}{1}.probe(1); % get the begin index of probe
                                    res_matrix(num_rec, 3) = res_cells{itest, iprobe, itime, igamma, idgamma, ilambda, idict}{1}.probe(prange); % get the end index of probe

                                    ptime = size(res_cells{itest, iprobe, itime, igamma, idgamma, ilambda, idict}{1}.time,2);
                                    res_matrix(num_rec, 4) = res_cells{itest, iprobe, itime, igamma, idgamma, ilambda, idict}{1}.time(1); % get the begin index of tiem
                                    res_matrix(num_rec, 5) = res_cells{itest, iprobe, itime, igamma, idgamma, ilambda, idict}{1}.time(ptime); % get the begin index of time

                                    res_matrix(num_rec, 6) = res_cells{itest, iprobe, itime, igamma, idgamma, ilambda, idict}{1}.gamma; % get gamma
                                    res_matrix(num_rec, 7) = res_cells{itest, iprobe, itime, igamma, idgamma, ilambda, idict}{1}.dictgamma; % get gamma
                                    res_matrix(num_rec, 8) = res_cells{itest, iprobe, itime, igamma, idgamma, ilambda, idict}{1}.lambda; % get lambda
                                    res_matrix(num_rec, 9) = size(res_cells{itest, iprobe, itime, igamma, idgamma, ilambda, idict}{1}.initDict, 2); % get dicctionary size

                                    res_matrix(num_rec, 10) = iiter; % get the inner iteration number

                                    res_matrix(num_rec, 11) = res_cells{itest, iprobe, itime, igamma, idgamma, ilambda, idict}{2}(1, iiter); % get step-1 obj
                                    res_matrix(num_rec, 12) = res_cells{itest, iprobe, itime, igamma, idgamma, ilambda, idict}{2}(2, iiter); % get step-2 obj
                                    res_matrix(num_rec, 13) = res_cells{itest, iprobe, itime, igamma, idgamma, ilambda, idict}{2}(3, iiter); % get train err
                                    res_matrix(num_rec, 14) = res_cells{itest, iprobe, itime, igamma, idgamma, ilambda, idict}{2}(4, iiter); % get test err

                                    num_rec = num_rec + 1;

                                end
                            end

                        end
                    end
                end
            end
        end
    end
end


    
end