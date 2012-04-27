function [ output_args ] = avg_plots( fs_avg, ksvd_avg, nmf_avg, filter )
%AVG_PLOTS Summary of this function goes here
%   Detailed explanation goes here
% plot the avg statistical result using bar chart
%   fs_avg, the feature sign average result
%   ksvd_avg, the ksvd result
%   nmf_avg, the non-negateive matrix factorization result


if filter
    res_filter = [1,1,1,0,0,0,1,1,1]';
    fs_avg(1:9,:) = fs_avg([1,4,5,6,7,8,9,2,3],:);
    ksvd_avg(1:9, :) = ksvd_avg([1,4,5,6,7,8,9,2,3],:);

    fs_avg_fil = fs_avg(logical(res_filter), :);
    ksvd_avg_fil = ksvd_avg(logical(res_filter), :);
    nmf_avg_fil = nmf_avg(logical(res_filter), :);
else
    fs_avg_fil = fs_avg;
    ksvd_avg_fil = ksvd_avg;
    nmf_avg_fil = nmf_avg;
    
end

% plot the accuracy
barinput = horzcat(fs_avg_fil(:,2), ksvd_avg_fil(:,2), nmf_avg_fil(:,2));
figure;
bar(cell2mat(barinput));
xlabel('dataset', 'fontsize', 12);
ylabel('accuracy', 'fontsize', 12);
%legend('2fs', 'ksvd', 'nmf');
legend('2fs', 'ksvd', 'nmf', 'location', 'northoutside', 'orientation', 'horizontal');
legend boxoff;


%plot precision
precision = horzcat(fs_avg_fil(:, 5), ksvd_avg_fil(:, 5), nmf_avg_fil(:, 5));
figure;
bar(cell2mat(precision));
xlabel('dataset', 'fontsize', 12);
ylabel('precision', 'fontsize', 12);
%legend('2fs', 'ksvd', 'nmf');
legend('2fs', 'ksvd', 'nmf', 'location', 'northoutside', 'orientation', 'horizontal');
legend boxoff;


%plot mcc
mcc = horzcat(fs_avg_fil(:, 6), ksvd_avg_fil(:, 6), nmf_avg_fil(:, 6));
figure;
bar(cell2mat(mcc));
xlabel('dataset', 'fontsize', 12);
ylabel('MCC', 'fontsize', 12);
%legend('2fs', 'ksvd', 'nmf');
legend('2fs', 'ksvd', 'nmf', 'location', 'northoutside', 'orientation', 'horizontal');
legend boxoff;


end

