function [res_matrix] = mul_res_stat(matdir, algo)
%usage, e.g. mul_res_stat('*.mat', '2fs')
%[res_matrix] = mul_res_stat(matdir, algo)

matfiles = dir(matdir);
first = 1;
%load all the mat files and process cell result
for i = 1:numel(matfiles)
	outstr = sprintf('%s\n', matfiles(i).name);
	disp(outstr);

	clearvars -except matfiles i first res_matrix algo;
	loadstrc = load(matfiles(i).name, 'objs');
	res_cells = loadstrc.objs;

	if strcmp(algo, '2fs')
		rmatrix = fs_res_stat(res_cells);

	elseif strcmp(algo, 'nmf')
		rmatrix = nmf_res_stat(res_cells);

	elseif strcmp(algo, 'ksvd')
		rmatrix = ksvd_res_stat(res_cells);

	end

	if first
		res_matrix = rmatrix;
		first = 0;
	else
		res_matrix = vertcat(res_matrix, rmatrix);
	end

end

end
