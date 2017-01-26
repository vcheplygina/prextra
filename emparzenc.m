%EMPARZENC EM-algorithm for semi-supervised learning by parzenc
%
%   [W,V] = EMPARZENC(A,B,N,FID)
%   W = A*EMPARZENC([],B,N,FID)
%
% INPUT
%   A       Labeled dataset used for training
%   B       Additional unlabeled dataset
%   N       Number of smoothing parameter steps (default 1)
%   FID     File ID to write progress to (default [], see PRPROGRESS)
%
% OUTPUT
%   W      Trained classifier, based on A and B
%   V      Trained classifier based on A only
%
% DESCRIPTION
% Using the EM algorithm the PARZENC classifier is used iteratively
% on the joint dataset [A;B]. In EM each step the labels of A are reset
% to their initial values. Initial labels in B are neglected. They
% are iteratively updated as soft labels obtained by classifying B
% by the actual W. The EM algorithm is run for a fixed smoothing
% parameter of PARZENC. This is repeated for smaller smoothing 
% parameters in N steps, using harmonic interpolation between HL and HU,
% in which HL is the smoothing parameter estimate obtained from PARZENML
% applied to A and HU the estimate obtained from PARZENML applied to B.
% For N = 1, the average of HL and HU is used.
%
% SEE ALSO
% DATASETS, MAPPINGS, EMCLUST, EMC, PARZENC, PARZENML, PRPROGRESS

% Copyright: R.P.W. Duin, r.p.w.duin@prtools.org
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function [w,v] = emc(a,b,n,fid)
	if nargin < 4, fid = []; end
	if nargin < 3, n = 1; end
	if nargin < 2, b = []; end
	if nargin < 1 | isempty(a)
		w = mapping(mfilename,'untrained',{b,classf,labtype,fid});
		w = setname(w,'EMParzen CLassifier');
		return
	end

	if size(a,2) ~= size(b,2)
		error('Datasets should have same number of features')
	end
	
	c = getsize(a,3);
	epsilon = 1e-6;
	nlab = getnlab(a);
	lablist = getlablist(a);
	a = setlabels(a,nlab); 
	a = setlabtype(a,'soft');
	%ws = scalem([+a; +b],'variance');
	ws = unitm;
	a = a*ws;
	b = b*ws;
	lab = zeros(size(b,1),c);
	hl = parzenml(a);
	b = dataset(+b);
	hu = parzenml(b);
	hl = max(hl,hu) * 1.05;
	hu = min(hl,hu) * 0.95;
	if n == 1
		h = (hl+hu)/2;
	else
		dh = (log(hl) - log(hu))/(n-1);
		h = exp([log(hl):-dh:log(hu)]);
	end
	c = a;
	
	first = 1;
	for j=1:length(h)
		hh = h(j);
		prprogress(fid,['\nem_classifier optimization, h = ' num2str(hh) '\n'])	
		change = 1;
		while change > epsilon
			w = parzenc(c,hh);
			if first, v = w; first = 0; end
			d = b*w;
			labb = d*classc;
			change = mean(mean((+(labb-lab)).^2));
			lab = labb;
			b = setlabtype(b,'soft',lab);
			c = [a; b];
			prprogress(fid,'  change = %d\n', change)
		end
	end
	
	J = getlabels(w);
	w = ws*setlabels(w,lablist(J,:));
	v = ws*setlabels(v,lablist(J,:));
	