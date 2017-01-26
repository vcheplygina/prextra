%TESTMC More classifier tests
%
% E = TESTMC(A,W,TYPE)
% E = TESTMC(A*W,TYPE)

function e = testmc(a,w,varargin)

prtrace(mfilename);

if nargin < 1
	error('Insufficient parameters given')
end

isdataset(a);

if nargin == 1
	d = a;
	type = {'conf'};
elseif ismapping(w)
	d = a*w;
	if nargin < 3
		type = {'conf'};
	else
		type = varargin;
	end
elseif isstr(w)
	d = a;
	type = {w};
elseif iscell(w)
	d = a;
	type = w;
else
	error('Illegal call')
end

d = d*maxc;
I = matchlablist(getlablist(d),getfeatlab(d));
J = I(d.nlab); % J(i) is the column of the true label of object i
e = zeros(1,length(type));
m = size(d,1);
for j=1:length(type)
	switch(type{j})
		case 'conf'
            d = +normm(d,1);
			dd = d([1:m]'+(J-1)*m);
			e(j) = 1-mean(dd);
		otherwise
			error('Illegal type')
	end
end
			
			