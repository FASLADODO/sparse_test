function [stat_matrix] = stat_compute(res_matrix)
%compute the Accuracy, AUC, ROC, Sensitivity, Precision, MCC value of $res_matrix
%$res_matrix, is in the format of {strkey, TP, TN, FP, FN}
%stat_matrix, is in the format of {strkey, accuracy, sensitivity(TPR), FPR, precision(ppv), MCC}


stat_matrix = cell(1, 6);

for idx = 1:size(res_matrix, 1)
	strkey = res_matrix{idx, 1};
	TP = res_matrix{idx, 2};
	TN = res_matrix{idx, 3};
	FP = res_matrix{idx, 4};
	FN = res_matrix{idx, 5};
	P = TP + FN;
	N = FP + TN;
	Pi = TP + FP;
	Ni = FN + TN;
	
	%accuracy, acc = (TP + TN)/(P + N)
	acc = (TP + TN)/(P + N);

	%sensitivity, TPR = TP/(FP+FN)
	TPR = TP/(TP + FN);

	%FPR, FPR = FP / (FP + TN)
	FPR = FP/(FP+TN);

	%precision, precision = TP / (TP + FP)
	precision = TP/(TP+FP);

	%MCC, MCC = (TP * TN - FP * FN) / sqrt(P*N*P'*N')
	mcc = (TP * TN - FP * FN) / sqrt(P*N*Pi*Ni);
	
	stat_matrix(idx, :) = {strkey, acc, TPR, FPR, precision, mcc};
	idx = idx + 1;
end

end
