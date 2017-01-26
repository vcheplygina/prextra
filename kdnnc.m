%KDNNC 1-Nearest Neighbor Classifier using KDTree
%
%   W = KNNC(A,TYPE)
%
% INPUT
%   A      Dataset
%   TYPE  'fast' (default) or 'good'
%
% OUTPUT
%   W      NN classifier using KDTree for finding nearest neighbors
%
% DESCRIPTION  
% This is a KD-Tree nearest neighbor implementation for PRTools using the
% KD-Tree package by Pramod Vemulapalli, see 
% http://www.mathworks.com/matlabcentral/fileexchange/26649-kdtree-implementation-in-matlab
% It should be in the Matlab search path.
%
% SEE ALSO
% MAPPINGS, DATASETS, KNNC

% Copyright: R.P.W. Duin, r.p.w.duin@prtools.org
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands


function out = kdnnc(a,par)

	prtrace(mfilename);
	
	kdtreecheck;
	name = 'KDTree NN Classifier';
	% No input data, return an untrained classifier.
	if nargin < 2, par = 'good'; end
	if (nargin == 0) | (isempty(a))
		out = mapping(mfilename,'untrained',par);
		out = setname(out,name);
	elseif isdataset(a) & isstr(par) % training
		if isempty(strmatch(par,char('fast','good')))
			error('KDTree type should be either ''fast'' or ''good''')
		end
		islabtype(a,'crisp');
		isvaldfile(a,1,2); % at least 1 object per class, 2 classes
		a = testdatasize(a);
		a = testdatasize(a,'objects');
		a = seldat(a);    % get labeled objects only
		[m,k,c] = getsize(a);
		data.kdtree = kd_buildtree(+a,0);
		data.nlab = getnlab(a);
		data.type = par;
		out = mapping(mfilename,'trained',data,getlablist(a),k,c);
		out = setname(out,name);
	elseif nargin == 2 & ismapping(par)  % execution
		kdtree = getdata(par,'kdtree');
		nlab = getdata(par,'nlab');
		if strcmp(getdata(par,'type'),'fast')
			kdtype = @kd_closestpointfast;
		else
			kdtype = @kd_closestpointgood;
		end
		a_data = +a;
		% PRTools want posteriors. In this case these are just 0/1 's
		out = zeros(size(a,1),size(getlab(par),1));
		m = size(a,1);
		s = sprintf('Testing %i objects: ',m);
		prwaitbar(m,s);
		for i=1:m;
			prwaitbar(m,i,[s int2str(i)]);
			j = kdtype(kdtree,a_data(i,:));
			out(i,nlab(j)) = 1;
		end
		prwaitbar(0);
		out = setdat(a,out);
	end
		
	