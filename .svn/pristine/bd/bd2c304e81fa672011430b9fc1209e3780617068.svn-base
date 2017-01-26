%    out = tree_eval(w,x)
%
function out = tree_eval(w,x)

n = size(x,1);
out = zeros(n,1);

for i=1:n

	v=w;
	% if the first split is already solving everything (1 obj. per class)
	if isa(v,'double')
		out(i,1) = v;
	end
	while (out(i,1)==0)
		if (x(i,v.bestf)<v.bestt)
			v = v.l;
		else
			v = v.r;
		end
		if isa(v,'double')
			out(i,1) = v;
		end
	end
end

