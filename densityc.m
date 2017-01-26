%    w = densityc(a,u)
%
% Train a density per class in dataset A. The density is defined in the
% cell array U, where each element U{i} is an untrained density mapping.
% When just a single untrained mapping U is defined, the same model will
% be used for all classes.
%
% For instance:
% >> a = gendatb;
% >> u = scalem([],'variance')*parzenm;
% >> w = densityc(a,u)
% 
function w = densityc(a,u)

if nargin<2
	u = gaussm;
end
if nargin<1 || isempty(a)
	w = mapping(mfilename,{u});
	w = setname(w,'Densitybased classifier');
	return
end

if ~istrained(u)
	[n,d,c]=getsize(a);
	if ~isa(u,'cell')
		% make sure we have a mapping per class
		u = repmat({u},c,1);
	else
		if length(u)~=c
			error('I need %s untrained mappings in the cell array.',c);
		end
	end
	% train them:
	dens = cell(c,1);
	for i=1:c
		dens{i} = seldat(a,i)*u{i};
	end
	% save the stuff:
	W.dens = dens;
	W.prior = getprior(a,0);
	w = mapping(mfilename,'trained',W,getlablist(a),d,c);
	w = setname(w,'Density based classifier');
else
	% extract data:
	W = getdata(u);
	n = size(a,1);
	c = length(W.dens);
	out = zeros(n,c);
	for i=1:c
		out(:,i) = +(a*W.dens{i});
	end
	w = setdat(a,out,u);
end
