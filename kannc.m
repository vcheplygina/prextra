%KANNC 1-Nearest Neighbor Classifier using ANN Matlab Wrapper
%
%   W = KANNC(A,K,OPTION)
%
% INPUT
%   A      Dataset
%   K      Number of desired neighbours
%   OPTION Options for ANNQUERY
%
% OUTPUT
%   W      NN classifier using the ANN Query package
%
% DESCRIPTION  
% This is the nearest neighbor implementation for PRTools using the ANN
% Matlab Wrapper package. It should be in the path. If needed download it
% http://webscripts.softpedia.com/scriptDownload/ANN-MATLAB-Wrapper-Download-33976.html
%
% SEE ALSO
% MAPPINGS, DATASETS, ANNQUERY, KNNC

% Copyright: R.P.W. Duin, r.p.w.duin@prtools.org
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands


function out = kannc(a,knn,opt)

	prtrace(mfilename);
	
	annquerycheck;
	name = 'KANNC';
	% No input data, return an untrained classifier.
	if nargin < 3, opt = []; end
	if nargin < 2, knn = 1; end
	if (nargin == 0) | (isempty(a))
		out = mapping(mfilename,'untrained',{knn,opt});
		out = setname(out,name);
	elseif isdataset(a) & ~ismapping(knn) % training
		islabtype(a,'crisp');
		isvaldfile(a,1,2); % at least 1 object per class, 2 classes
		a = testdatasize(a);
		a = testdatasize(a,'objects');
		a = seldat(a);    % get labeled objects only
		[m,k,c] = getsize(a);
		v.data = a;
		v.knn  = knn;
		v.opt  = opt;
		out = mapping(mfilename,'trained',v,getlablist(a),k,c);
		out = setname(out,name);
		out = setbatch(out,0);
	elseif nargin == 2 & ismapping(knn)  % execution
		% to avoid confusion, rename trained mapping
		w = knn;
		v = +w;   % get datafield
		[k,c] = size(w);
		nlab = getnlab(v.data);
		J = annquery(+v.data',(+a)',v.knn)';
		n = size(J,1); % no of test objects
		J = nlab(J);
		if size(J,2) == 1
			out = zeros(c,n);
			out([1:c:n*c]+J'-1) = ones(1,n);
			out = out';
		else
			out = hist(J',[1:c])';
		end
		% Use Bayes estimators for posteriors
		out = (out+ones(size(out,1),c))/(v.knn+c);
		out = setdat(a,out,w);
	else
		error('Illegal input')
	end
		
	