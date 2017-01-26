%FEATSEL2 Pairwise feature selection for classification
% 
% [W,R] = FEATSEL2(A,CRIT,K,T,FID)
%
% INPUT	
%   A    Training dataset
%   CRIT Name of the criterion or untrained mapping 
%        (default: 'NN', i.e. the 1-Nearest Neighbor error)
%   K    Number of features to select (default: K = 0, return optimal set)
%   T    Tuning dataset (optional)
%   FID  File ID to write progress to (default [], see PRPROGRESS)
%
% OUTPUT
%   W    Output feature selection mapping
%   R    Matrix with step-by-step results
%
% DESCRIPTION
%
% Best pairs of features are evaluated inidividually, i.e. independently
% from other pairs. 
% 
% SEE ALSO 
% MAPPINGS, DATASETS, FEATEVAL, FEATSELF, FEATSELLR, FEATSEL,
% FEATSELO, FEATSELB, FEATSELI, FEATSELP, FEATSELM, PRPROGRESS

% Copyright: A. Harol, R.P.W. Duin, r.p.w.duin@prtools.org
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function [w,r] = featself(a,crit,ksel,t,fid)
		 
  prtrace(mfilename);
	
  if (nargin < 2) | isempty(crit)
    prwarning(2,'no criterion specified, assuming NN');
    crit = 'NN';
  end
  if (nargin < 3) | isempty(ksel)
    ksel = 0;
  end
  if (nargin < 4)
    prwarning(3,'no tuning set supplied (risk of overfit)');
    t = [];
  end
	if (nargin < 5)
		fid = [];
	end
	
	if nargin == 0 | isempty(a)
		% Create an empty mapping:
		w = mapping(mfilename,{crit,ksel,t});
	else	
		[m,k] = size(a);
		for j1=1:k
			for j2=j1+1:k
				J(j1,j2) = feateval(a(:,[j1,j2]),crit,t);
			end
		end
		w = pairs_sort(J, ksel);
		w = featsel(k,w);
		prprogress(fid,'featself  finished\n')
	end
	w = setname(w,'Pairwise FeatSel');

return

function mn = pairs_sort(J, MAX_F)      
% Artsiom Harol
	if MAX_F == 0, MAX_F = size(J,1); end
  ind2=find(J(:));
  a=J(ind2);
  [b,ind_sort]=sort(-a);
  %ind_sort=(flipud(ind_sort'))';
  ind2=ind2(ind_sort);
  [m,n]=ind2sub(size(J),ind2);
  
  mn=[m(:)';n(:)'];
 % mn=[m(1:MAX_F)';n(1:MAX_F)'];
  [umn] = unique(mn);
  for k=1:length(umn)
    ind = find(umn(k) == mn);
    mn(ind(2:end)) = [];
  end  
  mn = mn(1:MAX_F)';
return;  
