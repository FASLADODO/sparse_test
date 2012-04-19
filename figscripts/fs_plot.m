function fs_plot(avg_matrix)
%fs_plot(avg_matrix), plots the average result of acc, sensitivity, tpr, fpr, mcc, precision etc.
%avg_matrix, is in the structure of {strkey, accuracy, sensitivity(TPR), FPR, precision(ppv), MCC}
%strkey, 'pb1=1,pb2=100,t1=1,t2=1500,ga=0.01,dga=0.01,lbd=0.01,dsize=100'

	% select the data to plot
	strkey = sprintf('pb1=%d,pb2=%d,t1=%d,t2=%d,ga=%g,dga=%g,lbd=%g,dsize=%d', probe_begin, probe_end, time_begin, time_end, gamma, dictgamma, lambda, dictsize);

	for idx = 1:size(avg_matrix, 1)%number of rows
		orikey = avg_matrix{idx, 1};

	end

	% fix pb1, pb2, t1, t2, ga, dga, get the dynamical impacts of lbd & dsize
	spattern = sprintf('');


	% fix pb1, pb2, t1, t2, lbd, dsize, get the dynamical impacts of ga & dga



function sel_paramter()
	

end

end
