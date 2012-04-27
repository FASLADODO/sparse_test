function [ nmatrix ] = nagios_trim( rmatrix )
%NAGIOS_TRIM Summary of this function goes here
%  trim nagios dataset
%  rmatrix, the input nagios result matrix
%  nmatrix, trimed new nagios result matrix

threshold = 0.2;
rtolnum = size(rmatrix, 2);
uplimit = floor(rtolnum * threshold);
ridxs = zeros(size(rmatrix, 1),1);

for rowidx = 1:size(rmatrix, 1)
   nempty = numel(find(rmatrix(rowidx, :)==-10)); % get the number of empty reslut entries
   if nempty > uplimit
      ridxs(rowidx) = 1;
   end
end

nmatrix = rmatrix;
nmatrix(logical(ridxs),:) = [];

nmatrix(nmatrix>1)=1;
nmatrix(nmatrix<1)=0;

end