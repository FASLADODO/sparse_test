function [Dout, Objout] = l1ls_featuresign2 (params)
% The feature-sign search algorithm
% L1-regularized least squares problem solver
%
% This code solves the following problem:
% 
%    minimize_d 0.5*||y - A*d||^2 + lambda*sum_j(S_j*||d - d_j||^2) + gamma*||d||_1
%
% params.X, the transposed cofficient vector from previous step (feature sign for learning confficients)
% params.Y, training set, transpose of original example singnals Y
% params.gamma, regularized parameter for l1 constraint
% params.Dinit, initial dictionary matrix's transpose
% params.S, similarity matrix
% params.lambda, regularized parameter for l2 constraint on dictionary loss
% function
% Dout, the output dictionary after this iteration
% Objout, output objective function value after this iteration

% initialize the parameters
A = params.X;
Y = params.Y;
gamma = params.dictgamma;
Xinit = params.Dinit;
S = params.S;
epsilon = params.epsilon;
%S = S.*(ones(size(S))-eye(size(S)));
lambda = params.lambda;

warning('off', 'MATLAB:divideByZero');

use_Xinit= false;
if exist('Xinit', 'var')
    use_Xinit= true;
end

% initialize and pre-computation
Dout= zeros(size(A,2), size(Y,2));
AtA = A'*A;
AtY = A'*Y;

% rankA = rank(AtA);
% constrain the max selected number of non-zeros
rankA = min(size(A,1)-5, size(A,2)-5);
%rankA = 3;

%%%%
% todo;  add another loop here to iteratively update the D's columns until
% they are converged
%%%%

next_iteration = 1;
percentage = 0.001; % percentage to control the number of violated columns in one iteration
violated_set = ones(1, size(Xinit, 2));
max_local_iter = 5;

while next_iteration
   outstr=sprintf('step2 feature sign iteration for new loop, atoms updated: %d', sum(violated_set));
   disp(outstr);
   
        for i=1:size(Y,2) % size(Y, 2) is equal to columns of D', i.e. this is an iteration over columns
            if violated_set(i)
                %print sth when 100 samples reached
                if mod(i, 100)==0, fprintf('.'); end %fprintf(1, 'l1ls_featuresign: %d/%d\r', i, size(Y,2)); end

                % use the init X values, i.e. use initial dictionary
                if use_Xinit
                    idx1 = find(Xinit(:,i)~=0); 
                    maxn = min(length(idx1), rankA);
                    xinit = zeros(size(Xinit(:,i)));
                    xmatrix = Xinit;
                    % from the non-zero cofficients, select maximum $maxn entries
                    xinit(idx1(1:maxn)) =  Xinit(idx1(1:maxn), i);
                    [Dout(:,i), fobj]= ls_featuresign_sub (A, Y(:,i), AtA, AtY(:, i), gamma, S, lambda, xinit, xmatrix, i);
                else
                    [Dout(:,i), fobj]= ls_featuresign_sub (A, Y(:,i), AtA, AtY(:, i), gamma, S, lambda);
                end
                %update dictionary in the same loop
                %Dout(:,i) = Dout(:,i)./norm(Dout(:,i));
                Xinit(:,i) = Dout(:,i);
                Objout = fobj;
            end
        end
   
    %%%
    %%% may put the normalization of dictionary here
    %%%
    % normalize dictionary after each iteration
    for ij = 1:size(Dout, 2)
        if norm(Dout(:,ij)) ~= 0 && violated_set(ij)
            Dout(:, ij) = Dout(:, ij)./norm(Dout(:,ij));
        end
    end
    
    
    violated_set = violated_matrix(xmatrix, Dout, epsilon); 
    proportion = sum(violated_set)/size(violated_set, 2);
    outstr = sprintf('proportion of atoms need update: %f, epsilon: %f, previous fobj: %f', proportion, epsilon, fobj);
    disp(outstr);
    max_local_iter = max_local_iter -1;
    if sum(violated_set) < percentage * size(violated_set, 2) || (max_local_iter <= 0) % boolean add?????
        next_iteration = 0;
    end
    
end

    
fprintf(1, '\n');

warning('on', 'MATLAB:divideByZero');

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% compute the columns' difference between two matrices
% two strategies could be used here
% 1. each time update all columns of D, until all columns are converged
% 2. each time only update those obviously different columns and their
% similar peer columns (from the similarity matrix), until no columns
% change obviously.
% we might need an epsilon to determin the convergence
function [violated_set] = violated_matrix(oldm, newm, epsilon)

diff = sqrt(sum((newm - oldm).^2, 1))./sqrt(sum(newm.^2, 1)); % check element-wise sqrt  sqrt(sum(S.^2,1));

%distance = 
outstr = sprintf('max diff: %f, min diff: %f', max(diff), min(diff));
disp(outstr);
violated_set = (diff>epsilon);

return;


function [x, fobj] = ls_featuresign_sub (A, y, AtA, Aty, gamma, S, lambda, xinit, xmatrix, cur_indx)
% xmatrix is Xinit, i.e. the initial dictionary

[L,M] = size(A); % return dim of A

rankA = min(size(A,1)-10, size(A,2)-10); % the max number of nonzero cofficients
%rankA =3;
%%%%%%
% Step 1: Initialize
%%%%%%
usexinit = false;
if ~exist('xinit', 'var') || isempty(xinit)
    xinit= [];
    x= sparse(zeros(M,1));
    theta= sparse(zeros(M,1));
    act= sparse(zeros(M,1));
    allowZero = false;
else
    % xinit = [];
    x= sparse(xinit);
    theta= sparse(sign(x)); %sign of the x_i
    act= sparse(abs(theta));% active set,
    usexinit = true;
    allowZero = true;
end

fname_debug = sprintf('../tmp/fsdebug_%x.mat', datestr(now, 30));

fobj = 0; %fobj_featuresign(x, A, y, AtA, Aty, gamma);


%%% optimality0 is condition(a), and optimality1 is condition(b)
ITERMAX=1000;
optimality1=false;
for iter=1:ITERMAX
    % check optimality0
    act_indx0 = find(act == 0); % get the zeros in active set
    
    %%%%%%%%
    % change the gradient below
    %%%%%%%%
    grad = AtA*sparse(x) - Aty + lambda.*matrxi_grad(xmatrix, cur_indx, S);
    %grad = AtA*sparse(x) - Aty;
    theta = sign(x);

    optimality0= false;
    
    %%%%%%%%
    % Step 2, from the zeros in active set, select the max gradient
    %%%%%%%%
    [mx,indx] = max(abs(grad(act_indx0)));

    % get max gradient which is greater than gamma into active set
    if ~isempty(mx) && (mx >= gamma) && (iter>1 || ~usexinit)
        act(act_indx0(indx)) = 1; % set the corresponding index to 1
        theta(act_indx0(indx)) = -sign(grad(act_indx0(indx))); % set the sign
        usexinit= false;
    else
        optimality0= true;
        if optimality1
	    %disp('line 194, break, optimality1 is true');
            break;
        end
    end
    act_indx1 = find(act == 1);  % get all the elements in active sets

    if length(act_indx1)>rankA
        warning('sparsity penalty is too small: too many coefficients are activated');
	%!!!!!! disp('sparsity penalty is too small: too many coefficients are activated');
        %return;
    end

    if isempty(act_indx1) %length(act_indx1)==0
        % if ~assert(max(abs(x))==0), save(fname_debug, 'A', 'y', 'gamma', 'xinit'); error('error'); end
        if allowZero, allowZero= false; continue, end
        return;
    end

    % if ~assert(length(act_indx1) == length(find(act==1))), save(fname_debug, 'A', 'y', 'gamma', 'xinit'); error('error'); end
    k=0;
    while 1
        k=k+1;

        if k>ITERMAX
            warning('Maximum number of iteration reached. The solution may not be optimal');
            %disp('Maximum number of iteration reached. The solution may not be optimal');
            % save(fname_debug, 'A', 'y', 'gamma', 'xinit');
            return;
        end

        if isempty(act_indx1) % length(act_indx1)==0
            % if ~assert(max(abs(x))==0), save(fname_debug, 'A', 'y', 'gamma', 'xinit'); error('error'); end
            if allowZero, allowZero= false; break, end
            return;
        end

        %%%%%%
        % Step 3: feature-sign step
        %%%%%%
        [x, theta, act, act_indx1, optimality1, lsearch, fobj] = compute_FS_step (x, A, y, AtA, Aty, theta, act, act_indx1, gamma, S, lambda, xmatrix, cur_indx);

        %%%%%%
        % Step 4: check optimality condition 1
        %%%%%%
        if optimality1 
	    %disp('line 239. break, optimality1 is true');
            break; 
        end;
        if lsearch >0 %condition (a) is not satisfied
            continue; 
        end;

    end
end

if iter >= ITERMAX
    warning('Maximum number of iteration reached. The solution may not be optimal');
    % save(fname_debug, 'A', 'y', 'gamma', 'xinit');
end

%%%%%%%
% might need to change the following grad
%%%%%%%
if 0  % check if optimality ???
    act_indx1 = find(act==1);
    grad = AtA*sparse(x) - Aty;
    norm(grad(act_indx1) + gamma.*sign(x(act_indx1)),'inf')
    find(abs(grad(setdiff(1:M, act_indx1)))>gamma)
end

[fobj,fgrad] = fobj_featuresign(x, A, y, AtA, Aty, gamma, lambda, xmatrix, cur_indx, S);

return;

    
% compute the gradient for s_ij*(||d-dj||^2)
function [grad] = matrxi_grad(matrix, indx, similarity)
    %grad = (sum(matrix,2) - matrix(:,indx)*(size(matrix,2)-1)).*similarity(:,indx);

    grad = zeros(size(matrix,1),1);
    for i = 1:size(matrix, 2)
        if i ~= indx
            grad = grad + similarity(indx, i)*(matrix(:, indx) - matrix(:, i));
        end
    end
return;


% compute sum_j(S_ij*D'_j)
function [sumv] = sumvect(matrix, indx, S)
    %sumv = zeros(size(matrix,1),1);
    sumv = matrix*S(:,indx);% - matrix(:,indx)*S(indx,indx);
    %for ii = 1:size(matrix,2)
    %    if indx ~= ii
    %        sumv = sumv + S(indx,ii)*matrix(:,ii);
    %    end
    %end
return;


% compute the loss between dictionary columns
function [loss] = dictloss(matrix, indx, cur_vect, S)
    
    loss = sum(sum((matrix - repmat(cur_vect, 1, size(matrix,2))).^2*S(:,indx)));

return;


% copmute the feature sign step
function [x, theta, act, act_indx1, optimality1, lsearch, fobj] = compute_FS_step (x, A, y, AtA, Aty, theta, act, act_indx1, gamma, S, lambda, xmatrix, cur_indx)

x2 = x(act_indx1);
% A2 = A(:, act_indx1);
AtA2 = AtA(act_indx1, act_indx1);
theta2 = theta(act_indx1);

% call matlab optimization solver..
% 'A \ B' left division of A into B, which is roughly the same as
% inverse(A)*B
 

%%%%%%%
% compute the x_new below
%%%%%%%
x_new = (AtA2 + lambda*(sum(S(cur_indx,:)))*speye(size(AtA2,1))) \ ( Aty(act_indx1) + lambda*sumvect(xmatrix(act_indx1,:), cur_indx, S) - gamma.*theta2 ) ;
% x_new = AtA2 \ ( Aty(act_indx1) - gamma.*theta2 ); % RR, get the accurate solution of quadratic problem
% opts.POSDEF=true; opts.SYM=true; % RR
% x_new = linsolve(AtA2, ( Aty(act_indx1) - gamma.*theta2 ), opts); % RR
optimality1= false;
if (sign(x_new) == sign(x2)) % signs are the same
    optimality1= true;
    x(act_indx1) = x_new;
    fobj = 0; %fobj_featuresign(x, A, y, AtA, Aty, gamma);
    lsearch = 1;
    return; 
end

% do line search: x -> x_new
progress = (0 - x2)./(x_new - x2);
lsearch=0;
%disp('start line search!');

%a= 0.5*sum((A2*(x_new- x2)).^2);

%%%%%%%
% !!! change the $a and $b expression below
%%%%%%%

a= 0.5*sum((A(:, act_indx1)*(x_new- x2)).^2);
b= (x2'*AtA2*(x_new- x2) - (x_new- x2)'*Aty(act_indx1));
fobj_lsearch = gamma*sum(abs(x2)) + dictloss(xmatrix(act_indx1, :), cur_indx, x2, S);


[sort_lsearch, idx_lsearch] = sort([progress',1]);
remove_idx=[];
for i = 1:length(sort_lsearch)
    t = sort_lsearch(i); 
    if t<=0 || t>1 
        continue; 
    end
    % get the step for next point
    s_temp= x2+ (x_new- x2).*t;
    
    %%%%%% 
    % compute the object function $fobj_temp below
    %%%%%%
    fobj_dict = lambda * dictloss(xmatrix(act_indx1,:), cur_indx, s_temp, S);
    fobj_temp = a*t^2 + b*t + gamma*sum(abs(s_temp)) + fobj_dict;
    
    %fobj_temp = a*t^2 + b*t + gamma*sum(abs(s_temp)); % get objective value at this point t 
    
    if fobj_temp < fobj_lsearch
        fobj_lsearch = fobj_temp;
        lsearch = t;
        if t<1  
            remove_idx = [remove_idx idx_lsearch(i)]; 
        end % remove_idx can be more than two..
    elseif fobj_temp > fobj_lsearch
        break;
    else
        if (sum(x2==0)) == 0
            lsearch = t;
            fobj_lsearch = fobj_temp;
            if t<1  
                remove_idx = [remove_idx idx_lsearch(i)]; 
            end % remove_idx can be more than two..
        end
    end
end

% if ~assert(lsearch >=0 && lsearch <=1), save(fname_debug, 'A', 'y', 'gamma', 'xinit'); error('error'); end

if lsearch >0
    % update x
    x_new = x2 + (x_new - x2).*lsearch;
    x(act_indx1) = x_new;
    theta(act_indx1) = sign(x_new);  % this is not clear...
end

% if x encounters zero along the line search, then remove it from
% active set
if lsearch<1 && lsearch>0
    %remove_idx = find(x(act_indx1)==0);
    remove_idx = find(abs(x(act_indx1)) < eps);
    x(act_indx1(remove_idx))=0;

    theta(act_indx1(remove_idx))=0;
    act(act_indx1(remove_idx))=0;
    act_indx1(remove_idx)=[];
end
%fobj_new = 0; %fobj_featuresign(x, A, y, AtA, Aty, gamma);
fobj = fobj_featuresign(x, A, y, AtA, Aty, gamma, lambda, xmatrix, cur_indx, S);
%fobj = fobj_new;

return;


% get the new objective function
function [func, gradient] = fobj_featuresign(x, A, y, AtA, Aty, gamma, lambda, xmatrix, cur_indx, S)

func= 0.5*norm(y-A*x)^2;
func= func+ gamma*norm(x,1);
func= func+ lambda*dictloss(xmatrix, cur_indx, x, S);

if nargout >1
    gradient= AtA*x - Aty;
    gradient= gradient+ gamma*sign(x);
end

outstr = sprintf('feature sign value is: %f', func);
%disp(outstr);

return;



function retval = assert(expr)
retval = true;
if ~expr 
    % error('Assertion failed');
    warning('Assertion failed');
    retval = false;
end

return

