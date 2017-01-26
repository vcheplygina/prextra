%
%    w = tree_train(x,y,opt)
%
function w = tree_train(x,y,opt)

% how good are we in this node?
err = tree_gini(y,opt.K);
if (err==0)

	w = y(1); % just predict this label

else
	% we split further
	n = size(x,1);

	% optionally, choose only from a subset
	if (opt.featsubset>0)
		fss = randperm(size(x,2));
		fss = fss(1:opt.featsubset);
	else
		fss = 1:size(x,2);
	end

	% check each feature separately:
	besterr = inf; bestf = []; bestt = []; bestj = []; bestI = [];
	for i=fss
		% sort the data along feature i:
		[xi,I] = sort(x(:,i)); yi = y(I);
		% run over all possible splits:
		for j=1:n-1
			% compute the gini
			err = j*tree_gini(yi(1:j),opt.K) + (n-j)*tree_gini(yi(j+1:n),opt.K);
			% and see if it is better than before.
			if (err<besterr)
				besterr = err;
				bestf = i;
				bestj = j;
				bestt = mean(xi(j:j+1));
				bestI = I;
			end
		end
	end

	% store
	w.bestf = bestf;
	w.bestt = bestt;
	%  now find the children:
	w.l = tree_train(x(bestI(1:bestj),:),y(bestI(1:bestj)),opt);
	w.r = tree_train(x(bestI(bestj+1:end),:),y(bestI(bestj+1:end)),opt);
end





