%GENARCHE Generate classifier archetypical data
%
%		A = GENARCHE(CLASSF,N,SEED,VERSION)
%
% INPUT
%  CLASSF : Classifer (untrained or name in string)
%  N      : Number of objects per class, default [50 50]
%  SEED   : Initial state for random generator, default as it is
%  VERSION: Version, default is most recent one
%
% OUTPUT
%  A      : Dataset

function a = genarche(classf,n,seed,version)

if nargin < 4, version = 0; end
if nargin < 3, seed = []; end
if nargin < 2, n = 50; end

if ~isempty(seed)
	rand('state',seed);
	randn('state',seed);
end

if ismapping(classf) & isuntrained(classf)
	classname = getname(classf);
elseif isstr(classf)
	classname = classf;
else
		error('Classifier should be given by string or by untrained mapping')
end

switch classname
	case {'qdc','QDA'}
		a = gendath([n n]);
		a = a * [4 1; -1 4] ./ sqrt(17);
	case {'udc','UDA'}
		a = gendath([n n]);
	case {'ldc','fisherc','LDA'}
		a = gendatd([n n]);
	case {'nmc','Nearest-Mean'}
		a = gendats([n n]);
	case {'treec','Dec-Tree'}
		m = genclass(n,[2/3 1/3]);
		a = [rand(n,1)*3+2 rand(n,1)*5]; 
		a = [a; rand(m(1),1)*2 rand(m(1),1)*5; rand(m(2),1)*5 rand(m(2),1)+5];
		a = dataset(a,genlab([n n]),'prior',0);
	case {'LM-Neural-Net'}
		m = genclass(n,[2/3 1/3]);
		a = [rand(n,1)*3+2 rand(n,1)*5]; 
		a = [a; rand(m(1),1)*2 rand(m(1),1)*5; rand(m(2),1)*5 rand(m(2),1)+5];
		a = a*[1 1; -1 1];
		a = dataset(a,genlab([n n]),'prior',0);
	case {'lmnc','neurc'}
		state1 = rand('state');
		state2 = randn('state');
		rand('state',0); 
		randn('state',0); 
		a = gendatb([200,200],0.5); 
		w = lmnc(a,5,40);              % find a feasible LM network with 5 hu
		rand('state',state1); 
		randn('state',state2); 
        b = gauss(4*n,[10 10],[9 0; 0 9]);
        b = exp(b/5);
        b = b - repmat(mean(b)+[1 1],4*n,1);
		%b = gendatb([2*n 2*n],1);    % generate a too large dataset
		d = +(b*w);                    % classify
		u = rand(4*n,1);               % reverse labels at random according posteriors
		J = find(u < min(+d,[],2));
		[dd,lab0] = max(+d,[],2);
        lab = lab0;
		lab(J) = 3-lab0(J);
		L1 = find(lab == 1);
		L2 = find(lab == 2);
		a = dataset([b(L1(1:n),:);b(L2(1:n),:)],genlab([n n]),'prior',0);
	case {'parzenc','Parzen'}
		a = gendatb([n n],2);
	case {'knnc','K-NN'}
		a = gendatb([n n],2);
		a = exp(a/7);
	case {'nnc','1nnc','spirals','1-NN'}
		x = rand(n,1)*4*pi+0.25*pi; 
		a = [x.*x.*sin(x),x.*x.*cos(x)];
		y = rand(n,1)*4*pi+0.25*pi; 
		b = [-y.*y.*sin(y),-y.*y.*cos(y)];
		a = dataset([a;b],genlab([n n]),'prior',0);
	case {'svc_nu'}
		% make a linearly separable problem, where one of the classes has
		% a terribly long moon-shape tail
		n = round(n/2);
		S = cat(3,5*eye(2),7*eye(2));
		extrablob = gauss([n n],[-35 35; 15 0],S);
		a = [gendatd([n n],2,3.5); extrablob];
    case {'asym-linear','asym-linear4'}
        %a = gendatd([2*n 2*n],2,3);
        b = gendats([n n]);
        b1 = seldat(b,1) + repmat([0 10],n,1);
        b2 = seldat(b,2) + repmat([0 -10],n,1);
        a = gendat([a; b1; b2],[n n]);
    case {'asym-linear2'}
				a1 = [rand(n,1) rand(n,1)*10];
				a2 = [rand(20*n,1)*0.1-1 rand(20*n,1)*10];
				b1 = [rand(n,1) rand(n,1)*10];
				b2 = [rand(2*n,1)*100+1 rand(2*n,1)*10];
				%b3 = [rand(18*n,1)+100 rand(18*n,1)*10];
				a = [a1;a2;b1;b2]*[1 1;-1 1];
				a = [a; +gauss(18*n,[40,80],(10*eye(2)))];
        a = dataset(a,genlab([21*n 21*n]));
        a = setprior(a,0);
				a = gendat(a,[n n]);
    case {'asym-linear3'}
				a1 = [rand(n,1)*10 rand(n,1)*100];
				b1 = [rand(n,1)*10+9 rand(n,1)*100+80];
				a1 = dataset([a1;b1]*[1 1; -1 1],genlab([n n]));
				b2 = gauss(n,[-100 150],10*eye(2));
				a2 = gauss(n,[-50 25],10*eye(2));
				a2 = dataset([a2;b2],genlab([n n]));
        a = setprior([a1;a2],0);
				a = gendat(a,[n n]);
    case {'SVM-1','asym-linear4'}
				a1 = [[rand(20*n,1)*5 rand(20*n,1)*100]; +gauss(2*n,[8,95],4*eye(2))];
				b1 = [[rand(20*n,1)*5+11 rand(20*n,1)*100]; +gauss(2*n,[8,5],4*eye(2))];
				a = [a1;b1]*[2 0; 0 1];
				%a1 = dataset([a1;b1]*[1 1; -1 1],genlab([n n]));
				a = dataset(a*[1 1; -1 1],genlab([22*n 22*n]));
				%b2 = gauss(n,[-20 80],10*eye(2));
				%a2 = gauss(n,[-90 80],10*eye(2));
				%a2 = dataset([a2;b2],genlab([n n]));
        a = setprior(a,0);
        %a = setprior([a1;a2],0);
				a = gendat(a,[n n]);
    case {'Logistic','asym-linear5'}
				a1 = [gauss(n,[5,50],[1 0; 0 100]); +gauss(n,[-3,25],4*eye(2))];
				b1 = [gauss(n,[7.5,50],[1 0; 0 100]); +gauss(n,[16,75],4*eye(2))];
				%a1 = dataset([a1;b1]*[1 1; -1 1],genlab([n n]));
				a = dataset([a1;b1]*[1 1; -1 1],genlab([2*n 2*n]));
				%b2 = gauss(n,[-20 80],10*eye(2));
				%a2 = gauss(n,[-90 80],10*eye(2));
				%a2 = dataset([a2;b2],genlab([n n]));
        a = setprior(a,0);
        %a = setprior([a1;a2],0);
				a = gendat(a,[n n]);
    case('line-plane')
        a = [rand(n,1) 0.25*((rand(n,1)*2).^2)];
        a = [a; rand(n,2)*[1 0; 0 0.01]+repmat([0 -0.01],n,1)]*[1 1; -1 1];
        a = dataset(a,genlab([n n]));
	case {'rbnc','RB-Neural-Net'}
		a = gendatb([n n]);
	case 'River'
		a = genriver([n n],1,0.3);
	case {'diamonds','Dis-Rep-LC'}
		a = rand(n,2);
		m = genclass(n,0.25*ones(1,4));
		d = (sqrt(2)-1)/2;
		b = [];
		b = [b; randxy(m(1),[-d -d],[0 1])];
		b = [b; randxy(m(2),[-d 1],[1 d+1])];
		b = [b; randxy(m(3),[1 0],[1+d 1+d])];
		b = [b; randxy(m(4),[0 -d],[1+d 0])];
		a = dataset([a;b],genlab([n n]),'prior',0);
		a = a*[1 1; 1 -1];
	case {'rbsvc','RB-SVM'}
		a = circ(n,1,0.5);
		m = genclass(n,[2/3 1/3]);
		b = circ(m(1),sqrt(1.5),1);
		c = circ(m(2),0.5);
		a = dataset([a;b;c],genlab([n n]),'prior',0);
	case 'circ2'
		a = circ(n,1,0.5);
		m = genclass(n,[2/3 1/3]);
		b = circ(m(1),sqrt(1.5),1);
		c = circ(m(2),0.5);
		a = dataset([a;b;c],genlab([n n]),'prior',0);
		a(:,2) = 3*a(:,2);
		a = a*[1 1; 1 -1];
	case 'chess4'
		ma = genclass(n,ones(1,8)/8);
		mb = genclass(n,ones(1,8)/8);
		a = randxy(ma(1),[0 0],[1 1]);
		a = [a; randxy(ma(2),[2 0],[3 1])];
		a = [a; randxy(ma(3),[1 1],[2 2])];
		a = [a; randxy(ma(4),[3 1],[4 2])];
		a = [a; randxy(ma(5),[0 2],[1 3])];
		a = [a; randxy(ma(6),[2 2],[3 3])];
		a = [a; randxy(ma(7),[1 3],[2 4])];
		a = [a; randxy(ma(8),[3 3],[4 4])];
		b = randxy(mb(1),[1 0],[2 1]);
		b = [b; randxy(mb(2),[3 0],[4 1])];
		b = [b; randxy(mb(3),[0 1],[1 2])];
		b = [b; randxy(mb(4),[2 1],[3 2])];
		b = [b; randxy(mb(5),[1 2],[2 3])];
		b = [b; randxy(mb(6),[3 2],[4 3])];
		b = [b; randxy(mb(7),[0 3],[1 4])];
		b = [b; randxy(mb(8),[2 3],[3 4])];
		a = dataset([a;b],genlab([n n]),'prior',0);
		a = a*[1 1; 1 -1];
	case 'chess41'
		ma = genclass(n,ones(1,8)/8);
		mb = genclass(n,ones(1,8)/8);
		a = randxy(ma(1),[0 0],[1 1]);
		a = [a; randxy(ma(2),[2 0],[3 1])];
		a = [a; randxy(ma(3),[1 1],[2 2])];
		a = [a; randxy(ma(4),[3 1],[4 2])];
		a = [a; randxy(ma(5),[0 2],[1 3])];
		a = [a; randxy(ma(6),[2 2],[3 3])];
		a = [a; randxy(ma(7),[1 3],[2 4])];
		a = [a; randxy(ma(8),[3 3],[4 4])];
		b = randxy(mb(1),[1 0],[2 1]);
		b = [b; randxy(mb(2),[3 0],[4 1])];
		b = [b; randxy(mb(3),[0 1],[1 2])];
		b = [b; randxy(mb(4),[2 1],[3 2])];
		b = [b; randxy(mb(5),[1 2],[2 3])];
		b = [b; randxy(mb(6),[3 2],[4 3])];
		b = [b; randxy(mb(7),[0 3],[1 4])];
		b = [b; randxy(mb(8),[2 3],[3 4])];
		a = dataset([a;b],genlab([n n]),'prior',0);
		a(:,2) = 3*a(:,2);
		a = a*[1 1; 1 -1];
	case {'Naive-Bayes'}
		a = [[randn(n,1)/6-0.5; randn(n,1)/6+0.5] 2*rand(2*n,1)-1];
		b = [2*rand(2*n,1)-1 [randn(n,1)/6-0.5; randn(n,1)/6+0.5]];
		a = dataset([a;b],genlab([2*n 2*n]));
		a = setprior(a,0);
		a = gendat(a,[n n]);
	otherwise
		error(sprintf('%s is not implemented',classname))
end

if ismapping(classf)
	a = setname(a,[getname(classf) '']);
else
	a = setname(a,classf);
end

function r = randxy(n,x,y)

if nargin < 1, n = 100; end
if nargin < 2, x = [0 0]; end
if nargin < 3, y = x + [1 1]; end

r1 = rand(n,2) .* repmat(y-x,n,1);		
r = r1 + repmat(x,n,1);	

function x = circ(n,r1,r2)

if nargin < 3, r2 = 0; end
if nargin < 2, r1 = 1; end

m = ceil(2*n*(r1*r1)/(r1*r1 - r2*r2));
x = rand(m,2) - repmat([0.5 0.5],m,1);
x = x*r1*2;
d = sqrt(sum(x.*x,2));
J = find(d < r1 & d > r2);
x = x(J(1:n),:);

function x = gausst(n,u,s,t)
x = randn(n,1);
x = t*x.*exp(abs(t)*x);
x = x - mean(x) + u;
x = s * x ./ std(x);

