function g = tree_gini(y,K)

out = zeros(1,K);
for k=1:K
	out(k) = mean(y==k);
end

g = out*(1-out)';

