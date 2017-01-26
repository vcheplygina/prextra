function [B, ll] = loglc2(A, W)
%LOGLC2 Logistic Linear Classifier
% 
%   W = LOGLC2(A, L)
% 
% INPUT
%   A   Dataset
%   L   Regularization parameter (L2)
%
% OUTPUT
%   W   Logistic linear classifier 
%
% DESCRIPTION  
% Computation of the linear classifier for the dataset A by maximizing the
% L2-regularized likelihood criterion using the logistic (sigmoid) function.
% The default value for L is 0.
% 
%
%  SEE ALSO 
%  MAPPINGS, DATASETS, LDC, FISHERC


    name = ['loglc2'];
    
    if ~exist('minFunc', 'file')
        error('LOGLC2 requires the minFunc optimizer. Please download it from www.di.ens.fr/~mschmidt/Software/minFunc.html and add it to the Matlab path.');
    end
    
    % Handle untrained calls like W = loglc2([]);
    if nargin == 0 || isempty(A)
        B = prmapping(mfilename);
        if nargin < 2
            W = 0;
        end
        
        name = ['loglc2_' num2str(W)];
    
        B = setname(B, name); 
        return;
        
    % Handle training on dataset A (use A * loglc2, A * loglc2([]), and loglc2(A))
    elseif (nargin == 1 && isdataset(A)) || (isdataset(A) && isa(W, 'double'))
        if nargin < 2
            W = 0;
        end
        name = ['loglc2_' num2str(W)];
    
        
        islabtype(A, 'crisp');
        isvaldfile(A, 1, 2);
        A = testdatasize(A, 'features');
        A = setprior(A, getprior(A)); 
        [~, k, c] = getsize(A);
        
        % Train the logistic regressor
        [data.E, data.E_bias] = train_logreg(+A', getnlab(A)', W);
        B = prmapping(mfilename, 'trained', data, getlablist(A), k, c);
        B = setname(B, name);
        
    % Handle evaluation of a trained LOGLC2 W for a dataset A 
    elseif (isdataset(A) && ismapping(W)) || (isa(A,'double') && ismapping(W))
        
        % Evaluate logistic classifier
        [~, test_post] = eval_logreg(+A', W.data.E, W.data.E_bias);
        A = prdataset(A); 
        B = setdata(A, test_post', getlabels(W));
        ll = [];
    
    % This should not happen
    else
        error('Illegal call');
    end
end

function [E, E_bias] = train_logreg(train_X, train_labels, lambda, E_init, E_bias_init)

    % Initialize solution
    if ~iscell(train_X)
        D = size(train_X, 1);
    else
        D = 0;
        for i=1:length(train_X)
            D = max(D, max(train_X{i}));
        end
    end
    [lablist, foo, train_labels] = unique(train_labels);
    K = length(lablist);
    if ~exist('E_init', 'var') || isempty(E_init)
        E = randn(D, K) * .0001;
    else
        E = E_init; clear E_init
    end
    if ~exist('E_bias_init', 'var') || isempty(E_bias_init)
        E_bias = zeros(1, K);
    else
        E_bias = E_bias_init; clear E_bias_init;
    end
    
    % Compute positive part of gradient
    pos_E = zeros(D, K);
    pos_E_bias = zeros(1, K);
    if ~iscell(train_X)
        for k=1:K
            pos_E(:,k) = sum(train_X(:,train_labels == k), 2);            
        end
    else
        for i=1:length(train_X)
            pos_E(train_X{i}, train_labels(i)) = pos_E(train_X{i}, train_labels(i)) + 1;            
        end
    end
    for k=1:K
        pos_E_bias(k) = sum(train_labels == k);
    end
    
    % Perform learning using L-BFGS
    x = [E(:); E_bias(:)];
    options.Method = 'lbfgs';
    %options.Display = 'on'; %DXD: nooooo!
    options.Display = 'off';
    options.TolFun = 1e-4;
    options.TolX = 1e-4;
    options.MaxIter = 5000;   
    if ~iscell(train_X)
        x = minFunc(@logreg_grad, x, options, train_X, train_labels, lambda, pos_E, pos_E_bias);
    else
        x = minFunc(@logreg_discrete_grad, x, options, train_X, train_labels, lambda, pos_E, pos_E_bias);
    end        
    
    % Decode solution
    E = reshape(x(1:D * K), [D K]);
    E_bias = reshape(x(D * K + 1:end), [1 K]);
end


function [est_labels, posterior] = eval_logreg(test_X, E, E_bias)

    % Perform labeling
    if ~iscell(test_X)
        log_Pyx = bsxfun(@plus, E' * test_X, E_bias');
    else
        log_Pyx = zeros(length(E_bias), length(test_X));
        for i=1:length(test_X)
            for j=1:length(test_X{i})
                log_Pyx(:,i) = log_Pyx(:,i) + sum(E(test_X{i}{j},:), 1)';
            end
        end
        log_Pyx = bsxfun(@plus, log_Pyx, E_bias');
    end
    [~, est_labels] = max(log_Pyx, [], 1);
    
    % Compute posterior
    if nargout > 1
        posterior = exp(bsxfun(@minus, log_Pyx, max(log_Pyx, [], 1)));
        posterior = bsxfun(@rdivide, posterior, sum(posterior, 1));
    end
end


function [C, dC] = logreg_grad(x, train_X, train_labels, lambda, pos_E, pos_E_bias)
%LOGREG_GRAD Gradient of L2-regularized logistic regressor
%
%   [C, dC] = logreg_grad(x, train_X, train_labels, lambda, pos_E, pos_E_bias)
%
% Gradient of L2-regularized logistic regressor.


    % Decode solution
    [D, N] = size(train_X);
    K = numel(x) / (D + 1);
    E = reshape(x(1:D * K), [D K]);
    E_bias = reshape(x(D * K + 1:end), [1 K]);

    % Compute p(y|x)
    gamma = bsxfun(@plus, E' * train_X, E_bias');
    gamma = exp(bsxfun(@minus, gamma, max(gamma, [], 1)));
    gamma = bsxfun(@rdivide, gamma, max(sum(gamma, 1), realmin));
    
    % Compute conditional log-likelihood
    C = 0;
    for n=1:N
        C = C - log(max(gamma(train_labels(n), n), realmin));
    end
    C = C + lambda .* sum(x .^ 2);
    
    % Only compute gradient when required
    if nargout > 1
    
        % Compute positive part of gradient
        if ~exist('pos_E', 'var') || isempty(pos_E)
            pos_E = zeros(D, K);
            for k=1:K
                pos_E(:,k) = sum(train_X(:,train_labels == k), 2);
            end
        end
        if ~exist('pos_E_bias', 'var') || isempty(pos_E_bias)
            pos_E_bias = zeros(1, K);
            for k=1:K        
                pos_E_bias(k) = sum(train_labels == k);
            end
        end

        % Compute negative part of gradient    
        neg_E = train_X * gamma';
        neg_E_bias = sum(gamma, 2)';
        
        % Compute gradient
        dC = -[pos_E(:) - neg_E(:); pos_E_bias(:) - neg_E_bias(:)] + 2 .* lambda .* x;
    end    
end
