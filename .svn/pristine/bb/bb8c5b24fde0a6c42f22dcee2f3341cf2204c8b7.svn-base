%AUC Area-Under-the-Curve error estimator
%
%   E = AUC(A,W)
%   E = AUC(A*W)
%   E = A*W*AUC
%
% Returns the ROC area_under_the_curve error 
% for the test set A on the classifier W.
%
% Two class problems only

% Copyright: R.P.W. Duin, r.p.w.duin@prtools.org
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function e = auc(a,w)

	if nargin == 0 | isempty(a)
		e = mapping('auc','fixed');
		return;
	end
	
	if nargin == 1, d = a*classc; end
	if nargin == 2, d = a*w*classc; end
	
	m = size(d,1);
	c = getsize(d,3);
	if c ~= 2
		error('Routine implemented for two-class problems only')
	end
	
	nlab = getnlab(d);
 	s = +d; s = sort(s(:));
	thr = [s' 1]; 
	
	e1 = []; e2 = [];

	% NLAB_OUT will be one where B is larger than THR.
	I = matchlablist(d.lablist,d.featlab); % Correct possible changes class orders
	
	nlab_out = (repmat(+d(:,I(1)),1,m*c+1) > repmat(thr,m,1));

	% ERRS will be 1 where the numeric label is unequal to NLAB_OUT
	% (i.e., where errors occur).
	errs = (repmat((nlab==1),1,m*c+1) ~= nlab_out);
			
	% Split the cases where ERRS = 1 into cases for CLAS (E1) and all 
	% other classes (E2).

	e1 = mean(errs(find(nlab==1),:),1);
	e2 = mean(errs(find(nlab~=1),:),1);
	
	J = [1:m*c];
	e = (e1(J+1)-e1(J))*(e2(J) + e2(J+1))'/2;