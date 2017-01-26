%WLOGLC2 Weighted Logistic Linear Classifier
% 
%   W = WLOGLC2(A, L, V)
% 
% INPUT
%   A   Dataset
%   L   Regularization parameter (L2)
%   V   Weights
%
% OUTPUT
%   W   Weighted Logistic linear classifier 
%
% DESCRIPTION  
% Computation of the linear classifier for the dataset A by maximizing the
% L2-regularized likelihood criterion using the logistic (sigmoid) function.
% The default value for L is 0.
% 
%
%  SEE ALSO 
%  MAPPINGS, DATASETS, LDC, FISHERC

function [b, ll] = wloglc2(a, lambda, v)

name = 'Weighted Logistic2';
if ~exist('minFunc', 'file')
   error('LOGLC2 requires the minFunc optimizer. Please download it from www.di.ens.fr/~mschmidt/Software/minFunc.html and add it to the Matlab path.');
end
   
if nargin<3
   v = [];
end
if nargin<2 || isempty(lambda)
   lambda = 0;
end
if nargin == 0 || isempty(a)
   b = prmapping(mfilename,{lambda,v}); 
   b = setname(b, name); 
   return;
end

if ~ismapping(lambda)
   % training
   islabtype(a, 'crisp');
   isvaldfile(a, 1, 2);
   a = testdatasize(a, 'features');
   a = setprior(a, getprior(a)); 
   [n, k, c] = getsize(a);

   % fix the weights:
   if isempty(n)
      v = ones(n,1);
   end
 
   % normalize     
   v = n*v(:)./sum(v);
 
        
   % Train the logistic regressor
   [data.E, data.E_bias] = train_logreg(+a', getnlab(a)', lambda, v);
   b = prmapping(mfilename, 'trained', data, getlablist(a), k, c);
   b = setdata(b, v, 'sampleweights');
   b = setname(b, name);
        
else
   % Evaluate logistic classifier
   W = getdata(lambda);
   [~, test_post] = eval_logreg(+a', W.E, W.E_bias);
   a = prdataset(a); 
   b = setdata(a, test_post', getlabels(lambda));
   ll = [];
end
    
end

function [E, E_bias] = train_logreg(train_X, train_labels, lambda, v, E_init, E_bias_init)

   % Initialize solution
   D = size(train_X, 1);
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
   VX = repmat(v',D,1).*train_X;
   for k=1:K
        pos_E(:,k) = sum(VX(:,train_labels == k), 2);            
   end
   for k=1:K
      I = find(train_labels==k);
      pos_E_bias(k) = sum(v(I));
   end
    
   % Perform learning using L-BFGS
   x = [E(:); E_bias(:)];
   options.Method = 'lbfgs';
   %options.Display = 'on'; %DXD: nooooo!
   options.Display = 'off';
   options.TolFun = 1e-4;
   options.TolX = 1e-4;
   options.MaxIter = 5000;   
   x = minFunc(@logreg_grad, x, options, train_X, train_labels, lambda, v, pos_E, pos_E_bias);
    
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


function [C, dC] = logreg_grad(x, train_X, train_labels, lambda, v, pos_E, pos_E_bias)
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
       C = C - v(n).*log(max(gamma(train_labels(n), n), realmin));
   end
   C = C + lambda .* sum(x .^ 2);
    
   % Only compute gradient when required
   if nargout > 1
   
      % Compute positive part of gradient
      if ~exist('pos_E', 'var') || isempty(pos_E)
         pos_E = zeros(D, K);
         VX = repmat(v',D,1).*train_X;
         for k=1:K
              pos_E(:,k) = sum(VX(:,train_labels == k), 2);            
         end
      end
      if ~exist('pos_E_bias', 'var') || isempty(pos_E_bias)
         pos_E_bias = zeros(1, K);
         for k=1:K        
            I = find(train_labels==k);
            pos_E_bias(k) = sum(v(I));
         end
      end

      % Compute negative part of gradient    
      neg_E = (repmat(v',D,1).*train_X) * gamma';
      neg_E_bias = gamma*v;
       
      % Compute gradient
      dC = -[pos_E(:) - neg_E(:); pos_E_bias(:) - neg_E_bias(:)] + 2 .* lambda .* x;
   end    
end
