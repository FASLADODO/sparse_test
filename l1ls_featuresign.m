function [Xout, fobj] = l1ls_featuresign (A, Y, gamma, Xinit)
% The feature-sign search algorithm
% L1-regularized least squares problem solver
%
% This code solves the following problem:
% 
%    minimize_s 0.5*||y - A*x||^2 + gamma*||x||_1
% 
% The detail of the algorithm is described in the following paper:
% 'Efficient Sparse Codig Algorithms', Honglak Lee, Alexis Battle, Rajat Raina, Andrew Y. Ng, 
% Advances in Neural Information Processing Systems (NIPS) 19, 2007
%
% Written by Honglak Lee <hllee@cs.stanford.edu>
% Copyright 2007 by Honglak Lee, Alexis Battle, Rajat Raina, and Andrew Y. Ng

warning('off', 'MATLAB:divideByZero');

use_Xinit= false;
if exist('Xinit', 'var')
    use_Xinit= true;
end

Xout= zeros(size(A,2), size(Y,2));
AtA = A'*A;
AtY = A'*Y;

% rankA = rank(AtA);
rankA = min(size(A,1)-10, size(A,2)-10);

for i=1:size(Y,2)
    %print sth when 100 samples reached
    if mod(i, 100)==0, fprintf('.'); end %fprintf(1, 'l1ls_featuresign: %d/%d\r', i, size(Y,2)); end
    
    % use the init X values, might be useful when dealing with similar rows
    % in D
    if use_Xinit
        idx1 = find(Xinit(:,i)~=0);
        maxn = min(length(idx1), rankA);
        xinit = zeros(size(Xinit(:,i)));
        xinit(idx1(1:maxn)) =  Xinit(idx1(1:maxn), i);
        [Xout(:,i), fobj]= ls_featuresign_sub (A, Y(:,i), AtA, AtY(:, i), gamma, xinit);
    else
        [Xout(:,i), fobj]= ls_featuresign_sub (A, Y(:,i), AtA, AtY(:, i), gamma);
    end
end
fprintf(1, '\n');

warning('on', 'MATLAB:divideByZero');

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [x, fobj] = ls_featuresign_sub (A, y, AtA, Aty, gamma, xinit)

[L,M] = size(A); % return dim of A

rankA = min(size(A,1)-10, size(A,2)-10); % the max number of nonzero cofficients
%rankA =2;
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
    %grad2 = AtA
    grad = AtA*sparse(x) - Aty;
    theta = sign(x);

    optimality0= false;
    
    %%%%%%%%
    % Step 2, select the max gradient
    %%%%%%%%
    [mx,indx] = max (abs(grad(act_indx0)));

    % get max gradient which is greater than gamma into active set
    if ~isempty(mx) && (mx >= gamma) && (iter>1 || ~usexinit)
        act(act_indx0(indx)) = 1; % set the corresponding index to 1
        theta(act_indx0(indx)) = -sign(grad(act_indx0(indx))); % set the sign
        usexinit= false;
    else
        optimality0= true;
        if optimality1
            break;
        end
    end
    act_indx1 = find(act == 1);  % get all the elements in active sets

    if length(act_indx1)>rankA
        warning('sparsity penalty is too small: too many coefficients are activated');
        return;
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
        [x, theta, act, act_indx1, optimality1, lsearch, fobj] = compute_FS_step (x, A, y, AtA, Aty, theta, act, act_indx1, gamma);

        %%%%%%
        % Step 4: check optimality condition 1
        %%%%%%
        if optimality1 
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

if 0  % check if optimality
    act_indx1 = find(act==1);
    grad = AtA*sparse(x) - Aty;
    norm(grad(act_indx1) + gamma.*sign(x(act_indx1)),'inf')
    find(abs(grad(setdiff(1:M, act_indx1)))>gamma)
end

[fobj, fgrd] = fobj_featuresign(x, A, y, AtA, Aty, gamma);

return;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x, theta, act, act_indx1, optimality1, lsearch, fobj] = compute_FS_step (x, A, y, AtA, Aty, theta, act, act_indx1, gamma)

x2 = x(act_indx1);
% A2 = A(:, act_indx1);
AtA2 = AtA(act_indx1, act_indx1);
theta2 = theta(act_indx1);

% call matlab optimization solver..
% 'A \ B' left division of A into B, which is roughly the same as
% inverse(A)*B

%%%%%%%
% change the x_new below
%%%%%%%

x_new = AtA2 \ ( Aty(act_indx1) - gamma.*theta2 ); % RR, get the accurate solution of quadratic problem
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
%a= 0.5*sum((A2*(x_new- x2)).^2);

%%%%%%%
% !!! change the $a and $b expression below
%%%%%%%

a= 0.5*sum((A(:, act_indx1)*(x_new- x2)).^2);
b= (x2'*AtA2*(x_new- x2) - (x_new- x2)'*Aty(act_indx1));
fobj_lsearch = gamma*sum(abs(x2));


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
    %!!! change the object function $fobj_temp below
    %%%%%%
    
    fobj_temp = a*t^2 + b*t + gamma*sum(abs(s_temp)); % get objective value at this point t 
    
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
fobj_new = 0; %fobj_featuresign(x, A, y, AtA, Aty, gamma);

fobj = fobj_new;

return;

%%%%%%%%%%%%
% get the new objective function, !!! should be replaced by a new one
%%%%%%%%%%%%

function [func, gradient] = fobj_featuresign(x, A, y, AtA, Aty, gamma)

func= 0.5*norm(y-A*x)^2;
func= func+ gamma*norm(x,1);

if nargout >1
    gradient= AtA*x - Aty;
    gradient= gradient+ gamma*sign(x);
end

return;

%%%%%%%%%%%%%%%%%%%%%

function retval = assert(expr)
retval = true;
if ~expr 
    % error('Assertion failed');
    warning ('Assertion failed');
    retval = false;
end
return
