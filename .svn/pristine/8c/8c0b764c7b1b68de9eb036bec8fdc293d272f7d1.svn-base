%
%          W = DECTREEC(A,N)
%
% Train a decision tree on A, using random feature subsets of size N.
% When N=0, no feature subsets are used.
% The decision tree training grows a full tree (no pruning), by
% splitting a single feature using one threshold. For the splitting the
% optimal GINI index is used.
function w = dectreec(a,featsubset)

if nargin<2
	featsubset = 0;
end
if nargin<1 || isempty(a)
	w = mapping(mfilename,{featsubset});
	w = setname(w,'Decision tree');
	return
end

if ~ismapping(featsubset)
	y = getnlab(a);
	opt.K = max(y);
	opt.featsubset = featsubset;
	if exist('decisiontree')==3
		v = decisiontree(+a,y,opt.K,featsubset);
	else
		v = tree_train(+a,y,opt);
	end
	w = mapping(mfilename,'trained',v,getlablist(a),size(a,2),opt.K);
else
	v = getdata(featsubset);
	n = size(a,1);
	if exist('decisiontree')==3

		if ~isa(v,'double')
			error('This tree should have been trained with the C-code');
		end
		out = decisiontree(v,+a);
	else
		if ~isa(v,'cell')
			error('This tree should have been trained with the Matlab code');
		end
		out = tree_eval(v,+a);
	end
	out = accumarray([(1:n)' out],ones(n,1));

	w = setdat(a,out,featsubset);
end
return
	
	
