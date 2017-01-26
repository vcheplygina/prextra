%LASSOC 
%
%      W = LASSOC(X, LAMBDA)
%
% Train the LASSO classifier on dataset X. LAMBDA is the regularization
% parameter.

% w = lassoc(x, lambda)
function w = lassoc(varargin)

argin = shiftargin(varargin,'scalar');
argin = setdefaults(argin,[],1);

if mapping_task(argin,'definition')

   w = define_mapping(argin,'untrained','LASSO classifier');

elseif mapping_task(argin,'training')

   [x,lambda] = deal(argin{:});
    % Unpack the dataset.
    islabtype(x,'crisp');
    %isvaldset(x,1,2); % at least 1 object per class, 2 classes
    [n,k,c] = getsize(x);

    % make sure a bias is added:
    x = [x ones(n,1)];

    if c ~= 2  % two-class classifier:
        error('Only a two-class classifier is implemented');
    end

    if exist('lasso')==3 % we have a own compiled mex code 
       beta=-lasso(+x,3-2*getnlab(x),lambda);
    else % hope that we have a modern Matlab with stats. toolbox:
       if ~exist('lasso')
          error('Cannot find the function lasso.m.');
       end
       beta=-lasso(+x,3-2*getnlab(x),'Lambda',lambda);
    end

    % now find out how sparse the result is:
    nr = sum(abs(beta)>1.0e-8);

    % and store the results:
    W.beta = beta; % the ultimate weights
    W.nr = nr;
    w = prmapping(mfilename,'trained',W,getlablist(x),k,c);
    w = setname(w,'LASSO classifier (l=%f)',lambda);

else
    [x,lambda] = deal(argin{1:2});
    % Evaluate the classifier on new data:
    W = getdata(lambda);
    n = size(x,1);

    % make sure a bias is added:
    x = [x ones(n,1)];
    %go:
    out = x*W.beta;

    % and put it nicely in a prtools dataset:
    w = setdat(x,sigm([-out out]),lambda);

end

return
