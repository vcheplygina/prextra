%GENRIVER Generate 2D 2-class River data
%
%		A = GENRIVER(N,L,S)
%
% INPUT
%   N    Number of objects in each of the classes (default: [50 50])
%   L    Length of the river, default: L = 1 (i.e. 2 bends)
%   S    Width of uniform class distribution on the river borders,
%        default: S = 1;
%
% OUTPUT
%   A    Generated dataset
%
% Generate a separable 'bending river'-like 2D dataset of two classes.
% Class priors are equal.
%
% EXAMPLE
%  repset = genriver([10 10],4);
%  trainset = genriver([],4);
%  w = proxm(repset,'d');
%  v = w*fisherc(trainset*w);
%  figure; scatterd(trainset); plotc(v)
%
% SEE ALSO
% DATASETS, PRDATASETS

% Copyright: A. Harol, R.P.W. Duin, r.p.w.duin@prtools.org
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function a = genriver(N,len,s,d)

	prtrace(mfilename);
	
	if nargin < 4, d = 0; end
	if nargin < 3, s = 1; end
	if nargin < 2 | isempty(len), len = 1; end
	if nargin < 1 | isempty(N), N = [50 50]; end

	p = [0.5 0.5];
	N = genclass(N,p);	
	x1=len*2*pi*rand(1,N(1));
	x2=sin(x1)+(s*rand(size(x1)))+d;
	d1 = prdataset([x1',x2'],ones(size(x1')));
	x3=len*2*pi*rand(1,N(2));
	x4 = sin(x3)+(-s+s*rand(size(x3)))-d;
	d2 = prdataset([x3',x4'],2*ones(size(x3')));
	a = setprior([d1;d2],p);
	
return
