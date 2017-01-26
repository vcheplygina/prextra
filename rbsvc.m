%RBSVC Automaric radial basis SVM
%
%   W = RBSVC(A,FID)

function w = rbsvc(a,fid)

if nargin < 2, fid = []; end
if nargin < 1 | isempty(a)
	w = mapping(mfilename,fid);
	return
end

nu = max(testk(a,1),0.05);

d = sqrt(+distm(a));
sigmax = min(max(d)); % max: smallest furthest neighboe distance

d = d + 1e100*eye(size(a,1));
sigmin = max(min(d)); % min: largest nearest neighbor distance

q = (1/sigmin - 1/sigmax)/9;
sig = [1/sigmax:q:1/sigmin];
sig = ones(1,length(sig))./[1/sigmax:q:1/sigmin];
if isempty(sig)
	sigopt = sigmin;
else
	w = [];
	errmin = inf;
	sigopt = 0;
	prprogress(fid,'     ')
	for j=1:length(sig)
		s = sprintf('\b\b\b\b\b%5.0f',j);
		prprogress(fid,s);
		err = crossval(a,svc_nu([],'r',sig(j),nu),10,1,fid);
		if err < errmin
			errmin = err;
			sigopt = sig(j);
		end
	end
	prprogress(fid,'\b\b\b\b\b')
end
w = svc_nu(a,'r',sigopt,nu);

