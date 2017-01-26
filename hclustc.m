%HCLUSTC Hierarchical clustering classifier
%
%    W = HCLUSTC(A,CTYPE,K)
%    W = HCLUSTC(A,CTYPE,K,U_D)
%
% It returns a hierarchical clustering classifier W, trained on dataset A.
% It is using CTYPE as the between-cluster distance or clustering
% criterion (see hclust.m for more details), and uses K clusters.
% This mapping generalizes to new, unseen data. This is done
% by computing the cluster distance between the new point and all the
% clusters obtained during training. The object is then assigned to the
% closest cluster.
% Per default the (squared) euclidean distance is used, but
% alternatively other proximity measures can be supplied in U_D.
%
% SEE ALSO
% HCLUST, PROXM

function w = hclustc(a,ctype,k,u_d)

if nargin<4
	u_d = proxm([],'d',2);
end
if nargin<3 || isempty(k)
	k = 10;
end
if nargin<2 || isempty(ctype)
	ctype = 's';
end
if nargin<1 || isempty(a)
	w = mapping(mfilename,{ctype,k,u_d});
	w = setname(w,'Hierarchical clustering (k=%d)',k);
	return
end

if ~ismapping(ctype)
	w = a*u_d;
	D = a*w;
	[lab, dendr] = hclust(D,ctype,k);
	x = dataset(+a,lab);

	W.u = u_d;
	W.x = x;
	W.k = k;
	W.ctype = ctype;
	w = mapping(mfilename,'trained',W,[],size(a,2),k);
	w = setname(w,'Hierarchical clustering (k=%d)',k);

else
	% evaluate the clustering on new data:
	W = getdata(ctype);
	n = size(a,1);
	out = zeros(n,W.k);
	% compute the distance to each cluster:
	for i=1:k
		w = seldat(W.x,i)*W.u;
		d = a*w;
		switch W.ctype
		case {'s','single'}
			out(:,i) = min(d,[],2);
		case {'c','complete'}
			out(:,i) = max(d,[],2);
		case {'a','average'}
			out(:,i) = mean(d,2);
		end
	end

	w = setdat(a,-out,ctype);

end
