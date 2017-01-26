%MODECLUST Clustering by mode-seeking in feature space
% 
% 	[LAB,J] = MODECLUST(A,K)
% 
% INPUT
%   A       Dataset
%   K       Number of neighbours to search for local mode (default: 10)
%
% OUTPUT
%   LAB     Cluster assignments, 1..K
%   J       Indices of modal samples
%
% DESCRIPTION
% A K-NN modeseeking method is used to assign objects to their nearest mode.
% The nearest neighbor implementation uses the ANN Matlab Wrapper package. 
% It should be in the path. If needed download it from
% http://webscripts.softpedia.com/scriptDownload/ANN-MATLAB-Wrapper-Download-33976.html
%
% LITERATURE
% Cheng, Y. "Mean shift, mode Seeking, and clustering", IEEE Transactions
% on Pattern Analysis and Machine Intelligence, vol. 17, no. 8, pp. 790-799,
% 1995.
% 
% SEE ALSO
% MAPPINGS, DATASETS, KMEANS, HCLUST, KCENTRES, PROXM

% Copyright: R.P.W. Duin, r.p.w.duin@prtools.org
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands


function [assign,J] = modeclust(a,k)

	prtrace(mfilename);
	
	% ANNQUERY needs to be in the path
	annquerycheck;

	% prepare data
	if (nargin < 2), k = 10; end
	m = size(a,1);
	b = (+a)';
	
	% Run the k-NN search, indices in J, distances in d
	[J,d]= annquery(b,b,k);

	% density estimate from distance to most remote neighbor
	f = 1./(max(d,[],1)+realmin)';
	
	% Find local indices of local density maxima in neighbourhood.
	[dummy,I] = max(reshape(f(J),size(J)));

	% Translate back to indices in all the data. N will contain the
	% dataset index of their initial mode estimate in the K-neighbourhood.
	N = J(I+[0:k:k*(m-1)]);

	% Climb the mode
	% Re-assign samples to the sample temporily mode is assigned to.
	% Iterate until assignments don't change anymore. Samples that then point 
	% to themselves are modes; all other samples point to the closest mode.

	M = N(N);
	while (any(M~=N))
		N = M; M = N(N);
	end

	% Use renumlab to obtain assignments 1, 2, ... and the list of unique
	% assignments (the modes).
	[assign,J] = renumlab(M');

return
