%AUCLPC Area Under the Curce Linear Programming Classifier
%
% [W,R] = AUCLPC(A,C,rtype, par, unitnorm, usematlab)
%
%   A  Dataset
%   C  Trade off parameter, default C = 10
%   W  Linear Classifier
%   R  Indices of support objects
%
% For other parameters see AUCLPM.
% The classifier optimizes the AUC for two-class problems using AUCLPM.
% The best threshold for the training set A is found using
% the class priors. Multi-class problems are solved using the
% one-against-rest scheme of MCLASSC (MODE = 'single').
%
% REFERENCE
% D.M.J. Tax and C.J. Veenman, Tuning the hyperparameter of an AUC-optimized
%    classifier, Proc. BNAIC 2005, 2005.
%
% SEE ALSO
% DATASETS, MAPPINGS, AUCLPM, MCLASSC

% Copyright: D.M.J. Tax, R.P.W. Duin, r.p.w.duin@prtools.org
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function [w,r] = auclpc(a,C, rtype, par, unitnorm, usematlab)

if (nargin < 6)
	usematlab = 0;
end
if (nargin < 5)
	unitnorm = 1;
end
if (nargin < 4)
	par = 0.25;
end
if (nargin < 3)
	rtype = 'subk';
end
if (nargin < 2)
	prwarning(3,'Lambda set to ten');
	C = 10; 
end

if nargin < 1 | isempty(a)
    w = mapping(mfilename,{C,rtype,par,unitnorm,usematlab});
    w = setname(w,'AUC_LP');
    return
end

[m,k,c] = getsize(a);
if c > 2
    %w = mclassc(a,feval(mfilename));
    w = mclassc(a,auclpc);
    w = setname(w,'AUC_LP');
    N = [];
    for j=1:c
        N = [N w.data{j}.data{1}.data{1}.data];
    end
    r = unique(N);
    return
end

v = auclpm(a,C,rtype, par, unitnorm, usematlab);
nlab = getnlab(a);
p = getprior(a);
d = +(a*v);
[dd,J] = sort(d);
n1 = sum(nlab==1);
n2 = sum(nlab==2);
%e = n1-n2;
%J1 = cumsum(nlab(J) == 1);
%J2 = n2-cumsum(nlab(J) == 2);
e1 = p(1)*cumsum(nlab(J) == 1)/n1 + p(2)*(1 - cumsum(nlab(J) == 2)/n2);
e2 = p(2)*cumsum(nlab(J) == 2)/n2 + p(1)*(1 - cumsum(nlab(J) == 1)/n1);

[min1,k1] = min(e1);
[min2,k2] = min(e2);
%if min1 < min2          % NO!! We should assume that AUCLPM gives a proper
%    if k1 == length(d)  % and well defined class order: 1 -> -, 2 > +
%        w0 = d(J(k1)) + realmin;
%    else
%        w0 = (d(J(k1))+d(J(k1+1)))/2;
%    end
%    w = v.data.u - v.data.v;
%else
    if k2 == length(d)
        w0 = d(J(k2)) + realmin;
    else
        w0 = (d(J(k2))+d(J(k2+1)))/2;
    end
    w = v.data.v - v.data.u;
%end
%J = find(w~=0);
prec = 1e-7;   % this is really choosing the significant guys
ss = sum(abs(w));
J  = find(abs(w) > ss*prec);

w = featsel(k,J)*affine(w(J),w0,a(:,J),getlablist(a),k,c); 
w = setout_conv(w,1);
w = setname(w,'AUC_LP');
r = J;
    
    
