function [ outD, outX, obj1 ] = sparse_learning( params )
%SPARSE_LEARNING Summary of this function goes here
%   Detailed explanation goes here
% params.train, training data with vars as rows and timestamp as columns
% params.dictsize, dictionary size
% params.lambda, regularized factor for l1ls_featuresign2, on dictionary
% loss function
% params.gamma, regularized factor for l1 norm on singal x
% params.initDict, initial dictionary for learning
% params.initX, initial cofficients?
% params.iternum, number of maximum iterations
%
% optional parameters
% params.test, test data set, same structure as params.train

global logfile;
logfile = fopen('log.txt', 'a');
log2file(logfile, datestr(now));
log2file(logfile, '%%%%%%%%%%%%%%%%%%%%');
log2file(logfile, 'params');
probe_range = [params.probe(1), params.probe(size(params.probe, 2))];
time_range = [params.time(1), params.time(size(params.time, 2))];
test_time = [params.test_time(1), params.test_time(size(params.test_time, 2))];
pstr = sprintf('dictsize: %d, gamma: %f, lambda: %f, iternum: %d, probe_range: %s, time_range: %s, test_time: %s', size(params.initDict,2), params.gamma, params.lambda, params.iternum, num2str(probe_range), num2str(time_range), num2str(test_time));
log2file(logfile, pstr);

Y = params.train;
D = params.initDict;
gamma = params.gamma;
lambda = params.lambda;
iternum = params.iternum;
test = params.test;

%optional parameters
display_figure = params.figure;

%temp output vars, 1st row: step1 feature sign obj; 2nd row: step2 feature
%sign obj; 3rd row: learned rmse err; 4th row: test rmse err;
obj1 = zeros(4,params.iternum);


%step 1, learn the similarity matrix, euclidean or hamming
simatrix = squareform(pdist(Y, 'hamming'));

simatrix = (ones(size(simatrix))-eye(size(simatrix))) - simatrix;


for iter = 1:iternum
    % step2, learn the cofficients, i.e. sparse coding
    disp('sparse coding...');
    log2file(logfile, 'sparse coding...');
    
    [xout, fobj] = l1ls_featuresign(D, Y, gamma);
    obj1(1, iter) = fobj;
    
    tparams.Y = Y';
    tparams.X = xout';
    tparams.gamma = gamma;
    tparams.Dinit = D';
    tparams.S = simatrix;
    tparams.lambda = lambda;
    tparams.epsilon = params.epsilon;
    
    outstr = sprintf('obj for step 1 feature sign: %f', fobj);
    disp(outstr);
    log2file(logfile, outstr);
    
    % step3, update the dictionary
    disp('dictionary updating...');
    log2file(logfile, 'dictionary updating...');
    
    [newD, obj] = l1ls_featuresign2(tparams);
    obj1(2, iter) = obj;
    
    %update the new dictionary and use it in the next iteration
    D = newD';
    outD = D;
    outX = xout;
    err = compute_err(outD, outX, Y);
    obj1(3, iter) = err;
    outstr = sprintf('iteration: %f, learned rmse on train_set is: %f, obj is: %f', iter, err, obj);
    log2file(logfile, outstr);
    disp(outstr);
    %rmse = err;
    
    [testX, ftobj, terr] = test_dict(outD, test, gamma);
    obj1(4, iter) = terr;
    
    % draw the figure of tmp results
    if display_figure
        figure(1);
        colormap = [1 1 0; 1 0 0; 1 1 1];
        pcolor(outD);
        title('dictionary');
        shading flat;
        
        figure(2);
        pcolor(outX);
        shading flat;
        title('cofficients');
        
        figure(3);
        pcolor(Y);
        shading flat;
        title('original data');

        figure(4);
        tmp = outD * outX;
        pcolor(tmp);
        shading flat;
        title('sparese representation - before normalized');
        
        figure(5);
        threshold = 0.2;
        tmp(tmp>=threshold)=1;
        tmp(tmp<threshold)=0;
        pcolor(tmp);
        shading flat;
        title('sparse representation - after normalized');
        errs = sqrt(sum(sum((tmp - Y).^2))/numel(Y));
        strout = sprintf('normalized rmse on the train set: %f', errs);
        disp(strout);
        log2file(logfile, strout);
        
    end
    
end


% save currrent workspace into file
outfile = sprintf('%s/wp.mat', params.outpath);
if exist(params.outpath, 'dir') || mkdir(params.outpath)
    save(outfile);
end

%plot the obj function during the iteration
figure(6);
plot(obj1(1:2, 1:size(obj1,2))');
dh = legend('step1_obj', 'step2_obj');
set(dh,'FontSize',12);
xlabel('iteration', 'fontsize', 12);
ylabel('value', 'fontsize', 12);
title('feature sign obj value', 'fontsize', 12);

%plot the rmse for train and test data
figure(7);
plot(obj1(3:4, 1:size(obj1,2))');
dh = legend('train_rmse', 'test_rmse');
set(dh,'FontSize',12);
xlabel('iteration', 'fontsize', 12);
ylabel('rmse', 'fontsize', 12);
title('RMSE value on train & test data', 'fontsize', 12);


fclose(logfile);

% return the learned cofficients, objective function value of step1 feature
% sign on test set, and the test error
%
function [testX, ftobj, err]=test_dict(dict_learned, test, gamma)
    disp('test the learned dictionary');
    [testX, ftobj] = l1ls_featuresign(dict_learned, test, gamma);
    err = compute_err(dict_learned, testX, test);
    outstt = sprintf('current test err(rmse): %f, ftobj: %f', err, ftobj);
    disp(outstt);
    log2file(logfile, outstt);
    
end

end




