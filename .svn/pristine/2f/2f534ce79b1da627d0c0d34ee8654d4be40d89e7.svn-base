function w = incsvc(a,ktype,par,C)
%INCSVC Incremental support vector classifier
%
%     W = INCSVC(A,KTYPE,PAR,C)
%
% INPUT
%   A      Dataset
%   KTYPE  Type of the kernel (optional; default: 'p')
%    PAR   Kernel parameter (optional; default: 1)
%    C     Regularization parameter (optional; default: 1)
%
%  OUTPUT
%    W     Mapping: Support Vector Classifier
%
%  DESCRIPTION
%  Optimizes a support vector classifier for the dataset A by an
%  incremental procedure to perform the quadratic programming. The
%  classifier can be of one of the types as defined by PROXM. Default is
%  linear (TYPE = 'p', PAR = 1). The kernel computation is done by
%  INCKERNEL which is more lightweight than PROXM.  C < 1 allows for
%  more class overlap. Default C = 1.
%
%  See also ADD_OBJ_CL, SVC, INCKERNEL

% Copyright: D.M.J. Tax, D.M.J.Tax@prtools.org
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

if nargin < 4 | isempty(C)
        C = 1;
        prwarning(3,'Regularization parameter C set to 1\n');
end
if nargin < 3 | isempty(par)
        par = 1;
        prwarning(3,'Kernel parameter par set to 1\n');
end
if nargin < 2 | isempty(ktype)
        ktype = 'p';
        prwarning(3,'Polynomial kernel type is used\n');
end
if nargin < 1 | isempty(a)
        w = mapping(mfilename,{ktype,par,C});
        w = setname(w,'Inc. Support Vector Classifier');
        return;
end

if ~isa(ktype,'mapping') %training
	%Basic sizes::
	[N,k,c] = getsize(a);
	%kernel definition:
	kernel = 'inckernel';
	kpar.type = ktype;
	kpar.s = par;
	%setting for displaying: (I know you can program it shorter, but this is
	%more clear):
	if N>=1000
		dodisplay = 1;
	else
		dodisplay = 0;
	end

	if c==2
		% randomize the data:
		%I = randperm(N); a = a(I,:);
		a = unrandomize(a);

		% settings for the program:
		global X_incremental;
		X_incremental = +a;
		y = 3 - 2*getnlab(a);

		% here we go:
		alf = zeros(N,1);   % weights
		grad = [];          % the gradient of all seen objects
		setR = [];          % 'rest' set
		Kr = [];            % kernel matrix of the rest objects(RxS)
		setD = [];
		setS = [];
		setE = [];
		Ke = [];
		Ks = 0;
		b = 0;
		R = inf;
		tol = 1e-8;

		% startup:
		for c=1:N
			if dodisplay & mod(c,100)==0
				fprintf('%d/%d ',c,N);
			end
			dd_message(6,'ADD obj %d\n',c);
			add_obj_cl;
%fprintf('%d [%f',c,alf(setS(1)));
%for jj=2:length(setS)
%    fprintf(' %f',alf(setS(jj)));
%end
%fprintf(']\n');
		end

		% make the classifier
		W.K = kernel;
		W.par = kpar;
		J = [setS;setE];
		W.J = J;
		W.sv = X_incremental(J,:);
		W.alf = y(J).*alf(J,:);
		W.b = b;
		w = mapping(mfilename,'trained',W,getlablist(a),k,2);
		w = setname(w,'Inc. Support Vector Classifier');

	else   % multi-class classifier:
		
		w = mclassc(a,mapping(mfilename,{ktype,par,C}));
		
	end

else   %execution
	W = +ktype;
	[n,d] = size(a);
	laba = getlab(a);
	orga = a;
	a = +a;
	global X_incremental;
	X_incremental = [W.sv; zeros(1,d)];
	nra = size(W.sv,1)+1; I = 1:nra;
	out = repmat(W.b,n,1);
	for i=1:n
		X_incremental(nra,:) = a(i,:);
		Ka = feval(W.K,W.par,nra,I);
		out(i) = out(i) + Ka(1:(nra-1))*W.alf;
	end
	newout = [out -out];
	w = setdat(orga,newout,ktype);
end

return



function a= unrandomize(a);
% Unrandomize a dataset;

[n,k,c]=getsize(a);
if c~=2
	error('I assume 2 classes');
end

nlab = getnlab(a);
I1 = find(nlab==1);
I2 = find(nlab==2);

if length(I1)<length(I2)
	cmin = length(I1);
else
	cmin = length(I2);
	tmpI = I1;
	I1 = I2;
	I2 = tmpI;
end

J=[];
J(1:2:2*cmin) = I1;
J(2:2:2*cmin) = I2(1:cmin);
J = [J,I2((cmin+1):end)'];
a = a(J,:);

return
