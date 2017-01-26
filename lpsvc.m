%LPSVC Sparse linear programming classifier according to Mangasarian
%
%   [W,E] = LPSVC(A,N,PoARAM)
%   W     = A*LPSVC([],N,PARAM)
%
% INPUT
%   A     Dataset
%   N     N-fold cross-validationocedure
%         N = 1: no cross-validation
%         N = size A: leave-one-out pr
%   PARAM Linear programming option
%         -1 Easy
%         0  Hard (default)
%         Otherwise: nu
%
% OUTPUT
%   W     Classifier
%   E     Internal error estimate
%
% DESCRIPTION
% This routine calls LPSVM by Fung and Mangasarian. If you use it, please
% refer to the below paper.
%
% REFERENCE
% G.M. Fung and O.L. Mangasarian, A Feature Selection Newton Method for
% Support Vector Machine Classification, Computational Optimization and
% Aplications, vol. 28, 2004, 185-202.
%
% SEE ALSO
% MAPPINGS, DATASETS, LINPROGC

% Elzbieta Pekalska, Robert P.W. Duin, e.pekalska@ewi.tudelft.nl
% Faculty of Electrical Engineering, Mathematics and Computer Science,
% Delft University of Technology, The Netherlands.

function [w,J] = lpsvm(a,n,param)
if nargin < 3, param =0; end
if nargin < 2, n =  1; end
name = 'LPSVC';
if nargin < 1 | isempty(a)
	w = mapping(mfilename,{n,param});
	w = setname(w,name);
	return
end
islabtype(a,'crisp');
isvaldset(a,1,2); % at least 1 object per class, 2 classes
[m,k,c] = getsize(a);
nlab = getnlab(a);	
	% The SVC is basically a 2-class classifier. More classes are
	% handled by mclassc.
if c == 2   % two-class classifier
	y = 3 - 2*nlab;
	[v,v0,ea,et] = lpsvm(+a,y,n,param);
	J = find(v~=0);
	if isempty(J)
		prwarning(1,'No features found. Fisher classifier is trained.')
		wf = featsel(k,1:k);
		w = wf*fisherc(a);
	else
		wf = featsel(k,J);
		flab = getfeatlab(a);
		w = wf*affine(v(J),-v0,flab(J,:),getlablist(a),length(J),c);
		w = cnormc(w,a);
	end
	
else
	
	[w,J] = mclassc(a,mapping(mfilename,{n,param}));
		
end







