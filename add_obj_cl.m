%Auxiliary script for incsvc.
%I especially did not make it a function to avoid a ton of arguments to
%pass. Maybe it gives a small speedup...

%fprintf('Add object number %d\n',c);
% currently seen objects plus new object
newsetD = (1:c)';

% compute the kernel matrix entry for the new object:
% (let us assume this is a row vector)
K = feval(kernel,kpar,c,newsetD);
Kc = (y(c)*y(newsetD)').*K;

%  % first see if the object is accepted or not
%  % (note that I could also have looked at the gradient, but then I
%  % don't have a very honest 'distance to the boundary')
%  this_label = K(c) - b;
%  if ~isempty(setS)
%    this_label = this_label - K(setS)*(y(setS).*alf(setS));
%  end
%  if ~isempty(setE)
%    this_label = this_label - K(setE)*(y(setE).*alf(setE));
%  end

% Compute the gradient for the new object:
grad(c,1) = y(c)*b - 1; %DXD here is a change for SVC
if ~isempty(setS)
	grad(c,1) = grad(c,1) + Kc(setS)*alf(setS);
end
if ~isempty(setE)
	grad(c,1) = grad(c,1) + Kc(setE)*alf(setE);
end

% Right, here we go, can be add object c?:
if grad(c)>0 % it is already classified ok,
	if length(setR)>0
		setR = [setR;c];   % put in the 'rest' set
		Kr = [Kr; Kc(setS)];  % update Kr
	else
		setR = c;
		Kr = Kc(setS);
	end

else         % we have to work

	done = 0;  % to check if we are done
	nrloops=0;
	while ~done

		% compute beta, not for the new object:
		beta = zeros(length(setS)+1,1);
		beta = -R*[y(c); Kc(setS)'];

		if any(beta>1e100) & ~isempty(setS)
			disp('Beta is too large: do you have identical objects with different labels??'); keyboard;
		end

		% compute gamma, also for the new object:
		if isempty(setS)
			% here something fishy is going on, when there is just a single
			% object added. Then we cannot freely move alf, and we can only
			% move b. In this case, b is moved until alf(c) enters setS
			% (that means, the gradient becomes 0).
			gamma = y(c)*y(newsetD);
			duptonow = y(c)*b;
		else
			gamma = Kc';
			if ~isempty(setE)
				gamma(setE) = gamma(setE) + [y(setE) Ke]*beta;
			end
			if ~isempty(setR)
				gamma(setR) = gamma(setR) + [y(setR) Kr]*beta;
			end
			gamma(c) = gamma(c) + [y(c) Kc(setS)]*beta;
			gamma(setS) = 0;
			duptonow = alf(c);
		end

      % now we have to see how large deltaAc can become...
      % (1) check the own upper bound:
      if isempty(setS)
        deltaAcisC = inf; %because we're moving b, and not alf!
      else
        deltaAcisC = C;
      end
      % (2) check if own gradient becomes zero:
      warning off;
      deltaGcis0 = duptonow-grad(c)./gamma(c);
      warning on;
      % object moves the wrong way:
      deltaGcis0(deltaGcis0<duptonow) = inf;
      % (3) check upper bounds of the SVs:
      deltaupperC = inf;
      if ~isempty(setS)
        warning off;
        deltaupperC = duptonow + (C-alf(setS))./beta(2:end);
        warning on;
        % negative changes do not count:
        deltaupperC(deltaupperC<duptonow) = inf;
        % object moves the wrong way (or not at all):
        deltaupperC(beta(2:end)<=0) = inf;
        [deltaupperC,nrS_up] = min(deltaupperC);
      end
      % (4) check lower bounds of the SVs:
      deltalowerC = inf;
      if ~isempty(setS)
        warning off;
        deltalowerC = duptonow + -alf(setS)./beta(2:end);
        warning on;
        % negative changes do not count:
        deltalowerC(deltalowerC<=duptonow) = inf;
        % object moves the wrong way (or not at all)
        deltalowerC(beta(2:end)>=0) = inf;
        [deltalowerC,nrS_low] = min(deltalowerC);
      end
      % (5) check E gradients to become 0:
      deltaGeis0 = inf;
      if ~isempty(setE)
        warning off; % divide by 0 is taken care of later...
        deltaGeis0 = duptonow -grad(setE)./gamma(setE);
        warning on;
        deltaGeis0(deltaGeis0<duptonow) = inf;
        deltaGeis0(gamma(setE)<=0) = inf;
        [deltaGeis0,nrE_0] = min(deltaGeis0);
      end
      % (6) check R gradients to become 0:
      deltaGris0 = inf;
      if ~isempty(setR)
        warning off; % divide by 0 is taken care of later...
        deltaGris0 = duptonow -grad(setR)./gamma(setR);
        warning on;
        %deltaGris0(deltaGris0<duptonow) = inf;
        deltaGris0(deltaGris0<=duptonow) = inf;
        deltaGris0(gamma(setR)>=0) = inf;
        [deltaGris0,nrG_0] = min(deltaGris0);
      end

      % which one is the most urgent one?
      deltas = [deltaAcisC; deltaGcis0; deltaupperC; deltalowerC;...
                deltaGeis0; deltaGris0];
      [maxdelta,situation] = min(deltas);
% fprintf('Situation %d (max_delta=%f)\n',situation,maxdelta);
%
      % update the parameters
      if isempty(setS) % then we only change b:
        %disp('setS is empty!');
        b = y(c)*maxdelta;
      else
        alf(c) = maxdelta;
        alf(setS) = alf(setS) + (maxdelta-duptonow)*beta(2:end);
        b = b + (maxdelta-duptonow)*beta(1);
%fprintf(' b: %f -> %f\n',bold,b);
%alf'
%fprintf('after: sum alf = %f\n',sum(alf));
      end
      grad = grad + (maxdelta-duptonow)*gamma;
      
      % do consistency check:
      I = find(alf<tol);
      if ~isempty(I)
        J = find(alf(I)<-tol);
%        if ~isempty(J)
%          disp('One of the alpha''s became < 0!');
%			 alf(I(J))
%          keyboard;
%		  end
        alf(I) = 0;
      end

      % update the sets:
      %fprintf('%d: situation %d  \n',c,situation);
      switch situation
      case 1   % object c goes to setE
        alf(c) = C; % just to be sure...
        if size(Ke,1)==0, Ke = []; end % make it really empty
        Ke = [Ke; Kc(setS)];
        setE = [setE; c];
        done = 1;
      case 2   % object c goes to setS
        Ks = [Ks [y(c); Kc(setS)']; y(c) Kc([setS; c])];
        Ke = [Ke Kc(setE)'];
        Kr = [Kr Kc(setR)'];
        % update R 
        if isempty(setS) % compute it directly (to avoid the inf's)...
          R = [-Kc(c) y(c); y(c) 0];
        else
          R = change_R(R,+c,beta,gamma(c));
        end
        setS = [setS;c];
        done = 1;
      case 3  % a support object hits upper bound
        j = setS(nrS_up);             % number of the object
        alf(j) = C;                   % just to be sure
        if size(Ke,1)==0, Ke=[]; end  % make it really really empty
        Ke = [Ke; Ks(nrS_up+1,2:end)];  % update Ke
        setE = [setE;j];              % add to setE
        Ks(nrS_up+1,:) = [];            % update all K's
        Ks(:,nrS_up+1) = [];
        Ke(:,nrS_up) = [];
        if ~isempty(Kr), Kr(:,nrS_up) = []; end
        setS(nrS_up) = [];            % remove from setS
        R = change_R(R,-nrS_up,beta,gamma(j));
      case 4  % a support object hits lower bound
        j = setS(nrS_low);             % number of the object
        alf(j) = 0;                    % just to be sure
        if size(Kr,1)==0, Kr = []; end % make really empty
        Kr = [Kr; Ks(nrS_low+1,2:end)];  % update Kr
        setR = [setR;j];               % add to setE
        Ks(nrS_low+1,:) = [];            % update all K's
        Ks(:,nrS_low+1) = [];
        if ~isempty(Ke), Ke(:,nrS_low) = []; end;
        if ~isempty(Kr), Kr(:,nrS_low) = []; end;
        setS(nrS_low) = [];            % remove from setS
        R = change_R(R,-nrS_low,beta,gamma(j));
      case 5  % an error becomes a support object
        j = setE(nrE_0);              % number of the object
        % adding to setS, means that all kernels have to be computed:
        K = feval(kernel,kpar,j,newsetD);
        Kj = (y(j)*y(newsetD)').*K;
        % to update R, we have to have the beta of object j:
        betaj = zeros(length(setS)+1,1);
        betaj = -R*[y(j); Kj(setS)'];
        Ks = [Ks; y(j) Kj(setS)];     % add row to Ks
        Kr = [Kr Kj(setR)'];          % update Kr
        Ke = [Ke Kj(setE)'];          % update Ke
        Ke(nrE_0,:) = [];
        setE(nrE_0) = [];             % update setE
        setS = [setS;j];              % add to setS
        Ks = [Ks [y(j); Kj(setS)']];   % and the extra column for Ks
        if length(betaj)==1 % compute it directly (to avoid the inf's)...
          R = [-Kj(j) y(j); y(j) 0];
        else
          % to update R, we also have to have the gamma of object j:
          gammaj = Ks(end,:)*[betaj;1] ;
          R = change_R(R,+j,betaj,gammaj);
        end
      case 6  % an other object becomes a support object
        j = setR(nrG_0);              % number of the object
        % adding to setS, means that all kernels have to be computed:
        K = feval(kernel,kpar,j,newsetD);
        Kj = (y(j)*y(newsetD)').*K;
        % to update R, we have to have the beta of object j:
        betaj = zeros(length(setS)+1,1);
        betaj = -R*[y(j); Kj(setS)'];
        Ks = [Ks; y(j) Kj(setS)];     % add row to Ks
        Ks = [Ks [y(j); Kj([setS;j])']]; % and the extra column for Ks
        Ke = [Ke Kj(setE)'];          % update Ke
        Kr = [Kr Kj(setR)'];          % update Kr
        Kr(nrG_0,:) = [];
        setS = [setS;j];              % add to setS
        setR(nrG_0) = [];             % update setR
        if length(betaj)==1 % compute it directly (to avoid the inf's)...
          R = [-Kj(j) y(j); y(j) 0];
        else
          % to update R, we also have to have the gamma of object j:
          gammaj = Ks(end,:)*[betaj;1] ;
          R = change_R(R,+j,betaj,gammaj);
        end
      end
      %fprintf('end situation\n');
      %alf,b
		nrloops = nrloops+1;
		%fprintf('nrloops = %d\n',nrloops);
		if nrloops>50, done=1; end
      %keyboard
    end % changing alfa's till stable solution is found
end % check for gradient<=0


% now we can also compute the R, compute the output for the x(setS):
% But fortunately, this value is
%   R^2 = offset + b;
% If we ignore the offset, we just have to use b!


