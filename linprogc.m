%LINPROGC Sparce linear programming classifier
%
%   W = LINPROGC(A,N,TYPE,PARAM)
%   W = A*LINPROGC([],N,TYPE,PARAM)
%
% INPUT
%   A     Dataset
%   N     Degree of classifier, 1,2,3,  ...
%   TYPE  Linear programming option, see below
%   PARAM Additional parameter
%
% OUTPUT
%   W     Classifier
%
% The parameter TYPE stands for the possible LP-machines:
% 'SIMPLE'   -  the most simple formulation; no sparse solution, hence 
%               all features are used.
% 'STANDARD' -  minimization of the training misclassification errors; no 
%               sparse solution, hence all features are used.
% 'C-SPARSE' -  a sparse solution (feature selection); a formulation similar 
%               to the LP_1 SVC, adjusted for the purpose of feature selection. 
%               PARAM is a tradeoff parameter, as in the traditional SVC.
%               Usually, PARAM = 1 or more. Default PARAM = 1.
% 'NU-SPARSE' - a sparse solution (feature selection); a formulation similar 
%               to the LP_1 SVC, based on the paper of Graepel, Herbrich, Smola etc
%               'Classification on proximity data with LP-machines' adjusted
%               for the purpose of feature selection. 
%               PARAM is a tradeoff parameter, usually PARAM = 0.05 or 0.1. 
%               It is an upper bound on the misclassfied training objects.      
%               So, for well separable problems, PARAM = 0.01 or PARAM = 0.05.
%               Default PARAM is the leave-one-out 1-NN error of A.
%
% DESCRIPTION
% A linear programming solution, similar to the LP_1 SVC is found for the 
% polynomial classifier (degree N, default N=1).
% Sparse solutions based on TYPE = C-SPARSE and TYPE = NU-SPARSE reduce the
% number of polynomial terms, thereby performing an automatic feature
% selection.
%
% This routine is effectively a wrapper around FSELLPC
%
% SEE ALSO
% MAPPINGS, DATASETS, POLYC, FSELLPC

% Elzbieta Pekalska, Robert P.W. Duin, e.pekalska@ewi.tudelft.nl
% Faculty of Electrical Engineering, Mathematics and Computer Science,
% Delft University of Technology, The Netherlands.

function w = linprogc(a,n,type,param)
if nargin < 4, param =[]; end
if nargin < 3, type = 'standard'; end
if nargin < 2, n =  1; end
name = [type '-LP'];
if nargin < 1 | isempty(a)
	w = mapping(mfilename,{n,type});
	w = setname(w,name);
	return
end
vs = scalem(a,'variance');
b = a*vs;
switch type
	case {'simple', 'Simple', 'SIMPLE'}
		type = 'Simple';
		param = 1;
	case {'standard', 'Standard', 'STANDARD'}
		type = 'Standard';
		param = 1;
	case {'c-sparse', 'C-Sparse', 'c-Sparse'}
		type = 'c-Sparse';
		if isempty(param), param = 1; end
	case {'nu-sparse','nu-Sparse','\nu-Sparse'}
		type = 'nu-Sparse';
		if isempty(param), param = testk(b,1); end
	otherwise
		error('Wrong Type')
end
vp = fsellpc([],1,type,param);
w = polyc(b,vp,n,1); 
w = vs*w;
name = [type '-LP'];
w = setname(w,name);
return
