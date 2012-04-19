%this script is used to test the sparse feature sign algorithm

load('nagios.mat');
%%read the test data file
if ~exist('nagios_trim','var')
    return;
    % dat = csvread('2006.csv', 0,1);
end
dat = nagios_trim;

traindim = 1500;
testdim = 500;
tpoints = (1:traindim:3001);
probe_range = 1:93:94;


%%{
dictsizes = [50];
cofs_gammas = [0.01];
dict_gammas = [0.1];
lambdas = [0.01];

objs = cell(2, size(probe_range, 2)-1, size(tpoints, 2)-2, size(cofs_gammas,2), size(dict_gammas,2), size(lambdas,2), size(dictsizes,2));


%train_probe = (1:100);
%test_probe = (1:100);
%%{
for itest = 1:1
    
for kk = 1:(size(probe_range, 2)-1)
    train_probe = (probe_range(kk):(probe_range(kk+1)-1));
    test_probe = train_probe;

    for jj = 1:(size(tpoints, 2)-2)
        train_time = (tpoints(jj):(tpoints(jj+1)-1));
        test_time = (tpoints(jj+1):(tpoints(jj+1)+testdim-1));

        % range from [500, 2000], [1500, 3000]

        %train_time = (1:1500);
        train = dat(train_probe, train_time);
      %  train(train==1)=0;
      %  train(train==-1)=1;

        % test set

        %test_time = (1501:2000);
        test = dat(test_probe, test_time);
     %   test(test==1)=0;
     %   test(test==-1)=1;

        for igamma=1:size(cofs_gammas, 2)
            for idictgamma = 1:size(dict_gammas, 2)
                for ilambda=1:size(lambdas, 2)
                    for idict=1:size(dictsizes, 2)
                        %clearvars -except dat train test dictsizes gammas lambdas igamma ilambda idict train_probe train_time test_probe test_time tpoints traindim testdim itest objs probe_range kk jj;

                        dictsize = dictsizes(idict);
                        randdict = train(:, rand_int(1,size(train,2),dictsize,1,1,0));
                        params.probe = train_probe;
                        params.time = train_time;
                        params.test_probe = test_probe;
                        params.test_time = test_time;

                        params.test = test;
                        params.train = train;

                        params.initDict = randdict;
                        params.gamma = cofs_gammas(igamma);
                        params.dictgamma = dict_gammas(idictgamma);
                        params.lambda = lambdas(ilambda);
                        params.iternum = 2;
                        params.figure = 0;
                        params.epsilon = 0.1;

                        probestr = num2str([test_probe(1), test_probe(size(test_probe,2))]);
                        train_timestr = num2str([train_time(1), train_time(size(train_time,2))]);
                        test_timestr = num2str([test_time(1), test_time(size(test_time,2))]);

                        outpath = sprintf('tround-%d_dict-%d_gamma-%g_dictgma-%g_lambda-%g_probe-[%s]_train-[%s]_test-[%s]', itest, dictsizes(idict), cofs_gammas(igamma), dict_gammas(idictgamma), lambdas(ilambda), probestr, train_timestr, test_timestr);
                        params.outpath = outpath;


                        [outD, outX, obj_func, outTestX] = sparse_learning(params);
                        params.testX = outTestX;

                        objs{itest, kk, jj, igamma, idictgamma, ilambda, idict} = {params, obj_func, outD, outX};

                        save('result.nagios.mat', 'objs');
                        %tsavefigures(outpath, 1, 0);

                    end
                end
            end
        end

    end
end
end

%}

%{
dictsize = 50;
% initial the dictionary
params.test = test;
params.train = train;
train_time = (1:1500);
train = dat(1:100, train_time);
train(train==1)=0;
train(train==-1)=1;
randdict = train(:,rand_int(1,size(train,2),dictsize,1,1,0));

params.train = train;
params.initDict = randdict;
params.gamma = 0.01;
params.lambda = 0.1;
params.iternum = 15;
params.figure = 1;
params.outpath = 'test';

[outD, outX, rmse]=sparse_learning(params);

tsavefigures(outpath, 1);


%}

%exit;
