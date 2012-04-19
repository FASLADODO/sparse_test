function [errs] = tolplot( stat_matrix, type )
%TOLPLOT Summary of this function goes here
%   Detailed explanation goes here
% plot the data from data_proc.m, statistic data(mean) about train & test err on
% different piece of dataset
% <data_start_time, data_end_time, gamma, dictgamma, lambda, dictsize, inner_index, train_err, test_err>
% type, type of err to plot, either 'train' or 'test'

% get the unique datasets
datasets = unique(stat_matrix(:,1:2),'rows');

% get the unique params
params = unique(stat_matrix(:, 3:6),'rows');
line_styles = auto_style(size(params, 1));

% initial the err matrix
errs = [];
 
num_figs = 1;

if strcmpi('train', type)
    % plot the figure for each dataset, test_err
    for idata =1:size(datasets, 1)
        start_idx = stat_matrix(:,1) == datasets(idata,1);
        end_idx = stat_matrix(:,2) == datasets(idata,2);
        data_idx = start_idx .* end_idx;
        % cur_dataset = stat_matrix(logical(data_idx), :);

        % for a dataset, generate a new figure
        figure(num_figs);
        num_figs = num_figs +1;
        train_err = [];
        legend_cell = {};

        for iparam = 1:size(params, 1)
            gamma_idx = stat_matrix(:, 3) == params(iparam, 1);
            dgamma_idx = stat_matrix(:, 4) == params(iparam, 2);
            lambda_idx = stat_matrix(:, 5) == params(iparam, 3);
            dict_idx = stat_matrix(:, 6) == params(iparam, 4);
            param_idx = gamma_idx .* dgamma_idx .* lambda_idx .* dict_idx;
            train_err_row = stat_matrix(logical(data_idx .* param_idx), 8); % get the train err
            train_err = [train_err, train_err_row];
            str_legend = sprintf('ga=%g,dga=%g,la=%g,di=%g', params(iparam,1), params(iparam,2), params(iparam,3), params(iparam, 4));
            legend_cell = horzcat(legend_cell, str_legend);
        end

        % plot(train_err, 'o-');
       
        hold on;
        for pidx = 1:size(train_err,2)
            plot(train_err(:,pidx), line_styles{pidx});
            disp(line_styles(pidx));
        end
        hold off;
        
        
        legend(legend_cell);
        title('train rmse');
        xlabel('iteration');
        ylabel('rmse value');
    end
    
elseif strcmpi('test', type)

    
    % plot the figure for each dataset, train err
    for idata =1:size(datasets, 1)
        start_idx = stat_matrix(:,1) == datasets(idata,1);
        end_idx = stat_matrix(:,2) == datasets(idata,2);
        data_idx = start_idx .* end_idx;
        % cur_dataset = stat_matrix(logical(data_idx), :);

        % for a dataset, generate a new figure
        figure(num_figs);
        num_figs = num_figs +1;
        test_err = [];
        legend_cell = {};

        for iparam = 1:size(params, 1)
            gamma_idx = stat_matrix(:, 3) == params(iparam, 1);
            dgamma_idx = stat_matrix(:, 4) == params(iparam, 2);
            lambda_idx = stat_matrix(:, 5) == params(iparam, 1);
            dict_idx = stat_matrix(:, 6) == params(iparam, 4);
            
            
           % gamma_idx = stat_matrix(:, 3) == 0.01;
           % dgamma_idx = stat_matrix(:, 4) == 0.01;
           % lambda_idx = stat_matrix(:, 5) == 0.1;
           % dict_idx = stat_matrix(:, 6) == 100;
            
            param_idx = gamma_idx .* dgamma_idx .* lambda_idx .* dict_idx;
          %   param_idx = lambda_idx;
            
            
            test_err_row = stat_matrix(logical(data_idx .* param_idx), 9); % get the test err
            test_err = [test_err, test_err_row];
             str_legend = sprintf('ga=%g,dga=%g,la=%g,di=%g', params(iparam,1), params(iparam,2), params(iparam,3), params(iparam, 4));
            % str_legend = sprintf('la=%g', params(iparam,1));
            
            legend_cell = horzcat(legend_cell, str_legend);
            
            %break;
        end

        %plot(test_err, '+-');
                
        % line_styles = auto_style(size(test_err, 2));
        % disp(size(test_err, 2));
        hold on;
        for pidx = 1:size(test_err,2)
            plot(test_err(:,pidx), line_styles{pidx});
            disp(line_styles(pidx));
        end
        hold off;
        %}
        
        %plot(test_err, 'r+-');
       
        
        legend(legend_cell);
        title('test rmse');
        xlabel('iteration');
        ylabel('rmse value');
    end
    
    
elseif strcmpi('err_sum', type)
    % sum the error on each dataset
    err_sum = [];
    
    for idata = 1:size(datasets, 1)
        start_idx = stat_matrix(:,1) == datasets(idata,1);
        end_idx = stat_matrix(:,2) == datasets(idata,2);
        data_idx = start_idx .* end_idx;
        % cur_dataset = stat_matrix(logical(data_idx), :);

        % for a dataset, generate a new figure
        test_err = [];
       
        %for iparam = 1:size(params, 1)
          %  gamma_idx = stat_matrix(:, 3) == params(iparam, 1);
          %  dgamma_idx = stat_matrix(:, 4) == params(iparam, 2);
          %  lambda_idx = stat_matrix(:, 5) == params(iparam, 1);
          %  dict_idx = stat_matrix(:, 6) == params(iparam, 4);
            
            gamma_idx = stat_matrix(:, 3) == 0.01;
            dgamma_idx = stat_matrix(:, 4) == 0.01;
            lambda_idx = stat_matrix(:, 5) == 0.1;
            dict_idx = stat_matrix(:, 6) == 100;
            iter_idx = stat_matrix(:, 7) == 10; % get the result of 10th iteration
            
            param_idx = gamma_idx .* dgamma_idx .* lambda_idx .* dict_idx .* iter_idx;
        
            test_err_row = stat_matrix(logical(data_idx .* param_idx), 9); % get the test err
            test_err = [test_err, test_err_row];
           % str_legend = sprintf('ga=%g,dga=%g,la=%g,di=%g', params(iparam,1), params(iparam,2), params(iparam,3), params(iparam, 4));
            % str_legend = sprintf('la=%g', params(iparam,1));
            
           % legend_cell = horzcat(legend_cell, str_legend);
            
            %break;
        %end
        err_sum = [err_sum, test_err];
        
    end
    
    errs = err_sum;
    
end


    
end

