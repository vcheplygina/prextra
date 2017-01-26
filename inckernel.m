function K = inckernel(par,I,J);
%INCKERNEL Kernel definition for incsvdd
%
%              K = INCKERNEL(PAR,I,J);
%
% Computation of the kernel function for the incremental SVDD. It is
% assumed that there is a global variable X_incremental, containing the
% objects. This is to avoid unnecessary overhead for very large
% datasets. Therefore we will also not use the 'standard' ways to comute
% the kernel (i.e. proxm).
%
% The kernel is defined by PAR.  We assume that PAR is a structure
% containing:
%  PAR.type   the kernel type
%       'kernel'       | 'k': X_incremental(I,J)
%       'polynomial'   | 'p': sign(xI*xJ'+1).*(xI*xJ'+1).^s
%       'exponential'  | 'e': exp(-(||xI-xJ||)/s)
%       'radial_basis' | 'r': exp(-(||xI-xJ||.^2)/(s*s))
%       'sigmoid'      | 's': sigm((sign(xI*xJ').*(xI*xJ'))/s)
%       'distance'     | 'd': ||xI-xJ||.^s
% And the special case:
%       'inline'       | 'i': use s as an inline function
%  PAR.s      the free parameter,
%             when more than 1 free parameter is required, then par.s
%             can be a vector of length 2 or more.
% The index vectors I and J indicate between which objects in
% X_incremental the kernel should be computed.
%
% See also: incsvdd, proxm

% Copyright: D. Tax, R.P.W. Duin, davidt@ph.tn.tudelft.nl
% Faculty of Applied Physics, Delft University of Technology
% P.O. Box 5046, 2600 GA Delft, The Netherlands

if nargin<2
	I = [];
	J = [];
end

global X_incremental;
if isempty(X_incremental)
  error('No data matrix X_incremental defined');
end
if isa(X_incremental,'dataset');
  error('Please make X_incremental a normal matlab array');
end

A = X_incremental(I,:);
B = X_incremental(J,:);

switch par.type
case {'kernel' 'k'}
	if isempty(I) & isempty(J)
		K = X_incremental;
	else
		K = X_incremental(I,J);
	end
case {'polynomial' 'p'}
	K = A*B'; 
	[n,d] = size(A);
	[m,d] = size(B);
	if par.s ~= round(par.s)
		K = K + ones(n,m);
		K = sign(K).*abs(K).^par.s;
	elseif par.s ~= 1
		K = K + ones(n,m);
		K = K.^par.s;
	end
case {'sigmoid' 's'}
	K = A*B'; 
	if length(par.s)>1
		K = sigm(K/par.s(1) + par.s(2));   %DXD: I need this sometimes
	else
		K = sigm(K/par.s);
	end
case {'exponential' 'e'}
	K = sqeucldistm(A,B);
	J = find(K<0);
	K(J) = zeros(size(J));
	K = exp(-sqrt(K)/par.s);
case {'radial_basis' 'r'}
	K = sqeucldistm(A,B);
	J = find(K<0);
	K(J) = zeros(size(J));
	K = exp(-K/(par.s*par.s));
case {'inv_radial_basis' 'i'}
	K = sqeucldistm(A,B);
	J = find(K<0);
	K(J) = zeros(size(J));
	K = exp(sqrt(K)/par.s);
case {'distance' 'd'}
	K = sqeucldistm(A,B);
	J = find(K<0);
	K(J) = zeros(size(J));
	if par.s ~= 2
		K = K.^(par.s/2);
	end
case {'inline' 'i'}
	error('Sorry, this does not work yet.');
case {'combp1s2'}
	K = sqeucldistm(A,B);
	J = find(K<0);
	K(J) = zeros(size(J));
	K = exp(-K/(2*2));
	K = par.s*A*B' + (1-par.s)*K; 
otherwise
	error(sprintf('Unknown proximity type: %s',par.type))
end

return
