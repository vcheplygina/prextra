function [B, ll] = randomforestc2(A, W)
%RANDOMFORESTC Random forest Classifier
% 
%   W = RANDOMFORESTC2(A, L)
% 
% INPUT
%   A   Dataset
%   L   Number of trees 
%
% OUTPUT
%   B   Random forest classifier 
%
% DESCRIPTION  
% 
%

    name = 'randomforestc2';
    if ~exist('classRF_train', 'file')
        error('Please download code from http://code.google.com/p/randomforest-matlab/');
    end
    
    % Handle untrained calls like W = randomforestc2([]);
    if nargin == 0 || isempty(A)
        if nargin < 2
            W = 100;
        end
        name=['randomforestc2_' num2str(W)];
     
        B = prmapping(mfilename); 
        B = setname(B, name); 
        return;
        
    % Handle training on dataset A (use A * loglc2, A * loglc2([]), and loglc2(A))
    elseif (nargin == 1 && isa(A, 'prdataset')) || (isa(A, 'prdataset') && isa(W, 'double'))
        if nargin < 2
            W = 100;
        end
        name=['randomforestc2_' num2str(W)];
     
        islabtype(A, 'crisp');
        isvaldfile(A, 1, 2);
        A = testdatasize(A, 'features');
        A = setprior(A, getprior(A)); 
        [~, k, c] = getsize(A);
        
        % Train the logistic regressor
        rf = classRF_train(+A, getnlab(A), W);
                
        B = prmapping(mfilename, 'trained', rf, getlablist(A), k, c);
        B = setname(B, name);
        
    % Handle evaluation of a trained RANDOMFORESTC2 W for a dataset A 
    elseif (isa(A, 'prdataset') && isa(W, 'prmapping')) || (isa(A, 'double') && isa(W, 'prmapping'))
        
        % Evaluate 
        [~, votes] = classRF_predict(+A, W.data);
           
        A = prdataset(A); 
        B = setdata(A, votes, getlabels(W));
        ll = [];
    
    % This should not happen
    else
        error('Illegal call');
    end
end
