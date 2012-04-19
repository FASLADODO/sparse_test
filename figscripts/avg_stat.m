function [avg_matrix] = avg_stat(stat_matrix)
%stat_matrix, is in the format of {strkey, accuracy, sensitivity(TPR), FPR, precision(ppv), MCC}
%avg_matrix, is in the format of {key, mAccuracy, mSensitivity, mFPR, mPrecision, mMCC};


keys = unique(stat_matrix(:,1));
avg_matrix = cell(1, 6);
indk = 1;
for ikey = 1:numel(keys)
	% get the mean vars below
	mAccuracy = mean(cell2mat(stat_matrix(strcmp(stat_matrix(:,1), keys(ikey)), 2)));
	mSensitivity = mean(cell2mat(stat_matrix(strcmp(stat_matrix(:,1), keys(ikey)), 3)));
	mFPR = mean(cell2mat(stat_matrix(strcmp(stat_matrix(:,1), keys(ikey)), 4)));
	mPrecision = mean(cell2mat(stat_matrix(strcmp(stat_matrix(:,1), keys(ikey)), 5)));
	mMCC = mean(cell2mat(stat_matrix(strcmp(stat_matrix(:,1), keys(ikey)), 6)));

	avg_matrix(indk, :) = {keys(ikey), mAccuracy, mSensitivity, mFPR, mPrecision, mMCC};
    indk = indk + 1;
    
end

end
