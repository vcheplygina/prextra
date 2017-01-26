function out = lessfx(par,x)
%LESSFX simple dataset mappings
%
%     PAR = LESSFX(TYPE,X)
%
% 'Train' or define a mapping of different types.
% TYPE can be:
%  1:    (x - m_2).^2 - (x - m_1).^2       normal nearest means
%  2:  ((x-m_2).^2)/s1 - ((x-m_1).^2)/s2   weighted nearest means
%  3:    (x-|m|_2).^2 - (x-|m|_1).^2       nearest medians
%  4: exp(-((x-M1).^2)./S1) - exp(-((x-M2).^2)./S2);
%  5:        "                 "     with medians instead of %  avergs
%
% New data is mapped using:
%     Y = LESSFX(PAR,X)
% where X in the input dataset, and PAR is obtained as above...
%
% This is used in LESS.M


if ~isstruct(par)  % we train the parameters of the function
	out = [];
	out.type = par;
	switch out.type
	case 0
		out.bla = [];
	case 1   % basic version
		out.u = +meancov(x);
	case 2   % mean-var
		[u,g] = meancov(x);
		out.u = +u;
		out.g(1,:) = diag(g(:,:,1)) + mean(diag(g(:,:,1)));
		out.g(2,:) = diag(g(:,:,2)) + mean(diag(g(:,:,2)));
	case 3   % median
		out.u(1,:) = med(+seldat(x,1));
		out.u(2,:) = med(+seldat(x,2));
	case 4   % mean-var
		[u,g] = meancov(x);
		out.u = +u;
		out.g(1,:) = diag(g(:,:,1)) + mean(diag(g(:,:,1)));
		out.g(2,:) = diag(g(:,:,2)) + mean(diag(g(:,:,2)));
	case 5   % median-MSD
		[u,g] = meancov(x);
		X1 = seldat(x,1);
		X2 = seldat(x,2);
		out.u(1,:) = med(+X1);
		out.u(2,:) = med(+X2);
		out.g(1,:) = medstd(+X1,out.u(1,:)) + mean(diag(g(:,:,1)));
		out.g(2,:) = medstd(+X2,out.u(2,:)) + mean(diag(g(:,:,2)));
	case 6 % median-MSD
		X1 = seldat(x,1);
		X2 = seldat(x,2);
		out.u(1,:) = med(+X1);
		out.u(2,:) = med(+X2);
		out.g(1,:) = medstd(+X1,out.u(1,:));
		out.g(2,:) = medstd(+X2,out.u(2,:));
	otherwise
		error('This function is not defined');
	end
else                 % we evaluate the function on new data:
	[m,k] = size(x);
	switch par.type
	case 0
		out = +x;
	case 1
		M1 = repmat(par.u(1,:),m,1);
		M2 = repmat(par.u(2,:),m,1);
		out = -(2*(+x).*(M2-M1) + M1.*M1 - M2.*M2);
	case 2
		M1 = repmat(par.u(1,:),m,1);
		M2 = repmat(par.u(2,:),m,1);
		S1 = repmat(par.g(1,:),m,1);
		S2 = repmat(par.g(2,:),m,1);
		out = ((x-M2).^2)./S2 - ((x-M1).^2)./S1;
	case 3
		M1 = repmat(par.u(1,:),m,1);
		M2 = repmat(par.u(2,:),m,1);
		out = -(2*(+x).*(M2-M1) + M1.*M1 - M2.*M2);
	case 4
		M1 = repmat(par.u(1,:),m,1);
		M2 = repmat(par.u(2,:),m,1);
		S1 = repmat(par.g(1,:),m,1);
		S2 = repmat(par.g(2,:),m,1);
		out = exp(-((x-M1).^2)./S1) - exp(-((x-M2).^2)./S2);
	case 5
		M1 = repmat(par.u(1,:),m,1);
		M2 = repmat(par.u(2,:),m,1);
		S1 = repmat(par.g(1,:),m,1);
		S2 = repmat(par.g(2,:),m,1);
		out = ((x-M2).^2)./S2 - ((x-M1).^2)./S1;
	case 6
		M1 = repmat(par.u(1,:),m,1);
		M2 = repmat(par.u(2,:),m,1);
		S1 = repmat(par.g(1,:),m,1);
		S2 = repmat(par.g(2,:),m,1);
		%out = ((x-M2).^2)./S2 - ((x-M1).^2)./S1;
		out = (sigm(10-((x-M1).^2)./(1*S1)) - sigm(10-((x-M2).^2)./(1*S2)));
	otherwise
		error('This function is not defined');
	end

end

return
