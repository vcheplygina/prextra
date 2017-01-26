function [lab, w] = mogm(a, k, share)
%MOGM Mixture-of-Gaussians clustering mapping
%
%   [LAB, W] = mogm(A, K, SHARE)
%
% INPUT
%   A       Dataset
%   K       Number of mixture components
%   SHARE   Use the same covariance for every component (true of false)
%
% OUTPUT
%   LAB     Cluster labels for training data
%   W       Clustering mapping


if isa(a,'double')
	a = prdataset(a,1);
end

    % Fitting of the model
    if isdataset(a) && ~ismapping(k)
        if nargin < 2
            error('Please input a dataset and the number of mixture components.');
        end
        if ~exist('share', 'var') || isempty(share)
            share = false;
        end
        tol = 1e-6;
        if share
            [lab, w] = emclust(a, ldc([], tol, tol), k, 'soft');
        else
            [lab, w] = emclust(a, qdc([], tol, tol), k, 'soft');
        end
        %[~, lab] = max(lab, [], 2);
        
    % Out-of-sample extension
    elseif isdataset(a) && ~ismapping(k)
        lab = a * k;
    else
        error('Erroneous call to MOGM.');
    end
