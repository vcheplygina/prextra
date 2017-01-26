function [w,I] = lessc(x, C, ftype, include_bias)
%LESSC Least Error in Sparse Subspaces classifier
%
%     W = LESSC(X, C, FTYPE, INCLUDE_BIAS)
%
% Train a linear classifier which also performs feature selection.
% In this version we do:
%                min \sum_i w_i + C*delta_i
%          s.t. forall_i   w^T f(x_i) > 1 - delta_i
%                       sum_i |w_i| = 1
% where f(x_i) is in principle free, but as a start we use the nearest
% mean idea:
%
%         f(x_i) = (x-mu2).^2 - (x-mu1).^2
% See for further definitions of f(x_i) lessfx.
%
% Dxd   15-3-2004
prtrace(mfilename);

if (nargin < 4)
	% To include a bias term in the model, we extend the number of features
	% by one:
	include_bias = 0;
end
if (nargin < 3)
	prwarning(3,'Use default function fx.');
	ftype = 1; 
end
if length(ftype)>1
	include_bias = ftype(2);
	ftype = ftype(1);
end
if (nargin < 2)
	prwarning(3,'C set to one');
	C = 1; 
end
if (nargin < 1) | (isempty(x))
	w = mapping(mfilename,{C,ftype,include_bias});
	w = setname(w,'LESS classifier');
	return
end

  	
if ~ismapping(C)   % train the mapping

	% Unpack the dataset.
	islabtype(x,'crisp');
	isvaldset(x,1,2); % at least 1 object per class, 2 classes
	[m,k,c] = getsize(x); 

	if c == 2   % two-class classifier

		% get -1/+1 labels:
		nlab = getnlab(x);
		y = 2*nlab-3;

		% train and apply the feature mapping:
		par = lessfx(ftype,x);
		f = lessfx(par,x);

		if (include_bias)
			f = [f ones(m,1)];
			k = k+1;
		end

		% In the LP formulation, we define the free parameter vector as:
		%  [delta; w]
		% setup the constraints:
		yf = -repmat(y,1,k).*f;
		% standard version when we have Ax<b  and  Aeq x = b;
		A = [-eye(m)    -(+yf)];
		b = -ones(m,1);
		%Aeq = [zeros(1,m) ones(1,k)]; beq = 1;
		%if (include_bias), Aeq(1,end)=0; end
		Aeq = []; beq = [];
		% function to optimize:
		c = [repmat(C,1,m) ones(1,k)];
		%c = [ones(1,m) repmat(C,1,k)];
		if (include_bias), c(end) = 0; end
		% upper and lower bounds:
		lb = zeros(m+k,1);
		if (include_bias), lb(end) = -inf; end
		ub = repmat(inf,m+k,1);

		% optimize
		if (exist('glpkmex')==3)
			[out,dummy]=glpkmex(1,c',A,b,repmat('U',m,1),lb,[],repmat('C',m+k,1));
		else
			out = linprog(c,A,b,Aeq,beq,lb,ub);
		end
		w = out((m+1):end);

		% find out how many features are relevant:
		if (include_bias)
			I = find(abs(w(1:(end-1)))>1e-8);
			nr = length(I);
		else
			I = find(abs(w)>1e-8);
			nr = length(I);
		end

		% Store the classifier
		W.extend = include_bias;
		W.par = par;
		W.w = w;
		W.nr = nr;
		w = mapping(mfilename,'trained',W,getlablist(x),size(x,2),2);
		w = setname(w,'LESS classifier');
	
	else   % multi-class classifier:
		
		%error('Multiclass not implemented yet');
		w = mclassc(x,mapping(mfilename,{C,ftype,include_bias}));
		v = w.data{1}.data{1}.data.w;
		for i=2:length(w.data)
			v = v + w.data{i}.data{1}.data.w;
		end
		I = find(abs(v)>0);
		
	end	
else
	% Evaluate the classifier on new data:
	W = getdata(C);

	% It is a simple linear classifier:
	if (W.extend)
		out = [lessfx(W.par,x) ones(size(x,1),1)]*W.w;
	else
		out = lessfx(W.par,x)*W.w;
	end

	% and put it nicely in a prtools dataset:
	w = setdat(x,sigm([out -out]),C);

end
		
return
