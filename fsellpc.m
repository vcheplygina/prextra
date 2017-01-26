%FSELLPC LP (linear programming) classifier
%
%   [W1,W2,W3] = fsellpc(A,BIAS,TYPE,PARAM)
%
% INPUT
%   A     Dataset
%   BIAS  YES or NO (optional; default: 1 (YES))
%   TYPE  Type of a classifier
%         'SIMPLE'   - the most simple formulation; all features used; PARAM = [];
%         'STANDARD' - minimization of the training misclassification errors;
%                      all features used; PARAM = [];
%         'C-SPARSE' - feature selection (sparse solution); a formulation similar
%                      to the LP_1 SVM; PARAM is a tradeoff parameter.
%                      (optional; DEFAULT: 1).
%         'MU-SPARSE' - feature selection (sparse solution); a formulation similar
%                       to the LP_1 SVM, based on the paper of Graepel, Herbrich, Smola etc
%                       'Classification on proximity data with LP-machines'.
%                       PARAM is a tradeoff parameter, usually PARAM = 0.05 or 0.1.
%                       It is an upper bound on the misclassfied training objects.
%                       So, for well separable problems, PARAM = 0.01 or PARAM = 0.02.
%                       (optional; DEFAULT: the LOO 1-NN error * 1.3).
%  PARAM Parameter connected to the TYPE, as above
%
% OUTPUT
%   W1  LP Classifier in the complete space
%   W2  LP Classifier in a reduced space
%   W3  Feature selection mapping; the indices of chosen features are in +W3.
%
% DEFAULTS
%   BIAS  = 1
%   TYPE  = 'STANDARD'
%   PARAM = []
%
% DESCRIPTION
% Classification problem on a N x M data A with LP-machines. Assume a two-class
% problem. Let DLPC select J features. Then:
% W1 is an M x 2 classifier in the original space, W2 is an J x 2 classifier
% in the feature space defined by the J chosen features and W3 is an M x R feature
% selection such that W1 = W3 * W2. Note that the indices of the selected features
% can be retrieved by +W3.
%
% A linear classifier is built on A:
%
%     f(A(x,F)) = A(x,F) * w + w0,
%
% where F are the features and w are the weights. If BIAS is 1, then w0 is
% also sought, otherwise it equals 0. This means that the hyperplane is
% forced to go through the origin.
%
% For C-class problems, C classifiers are trained, one against all others.
% In such a case, only W1 is returned and W3 in now NOT a feature selection,
% but directly the indices of the selected features.
%
% DEFAULT:
%   BIAS  = 1
%   TYPE  = 'STANDARD'
%   PARAM = 1
%

% Elzbieta Pekalska, Robert P.W. Duin, e.pekalska@tudelft.nl
% Faculty of Electrical Engineering, Mathematics and Computer Science,
% Delft University of Technology, The Netherlands.




function [W1,W2,W3] = fsellpc(a,is_w0,type,par,usematlab,prec)

if nargin < 6, prec = 1e-7; end
if nargin < 5, usematlab = 0; end
if nargin < 3 | isempty(type), type = 'standard'; end
if nargin < 4 | isempty(par),
  switch upper(type)
    case {'MU-SPARSE'}
      par = max(0.01,1.3*testk(a,1));     % upperbound error: 1.3 * loo 1-nn error
    case {'C-SPARSE'}
      par = 1;
    case {'SIMPLE','STANDARD'},
      par = [];
    otherwise
      disp(type)
      error('Wrong type.')
  end
end
if nargin < 2 | isempty(is_w0), is_w0 = 1; end
if nargin < 1 | isempty(a)
  W1 = mapping(mfilename,{is_w0,type,par,usematlab});
  W1 = setname(W1,'FSELLPC');
  W2 = [];
  W3 = [];
  return
end



if ~isdataset(a),
  error('The first parameter should be a dataset.')
end
if ~isnumeric(is_w0) | (is_w0 ~= 0 & is_w0 ~= 1),
  error('The second parameter should be 0 or 1.');
end


lab      = getnlab(a);
lablist  = getlablist(a);
[m,k,C]  = getsize(a);


z = (is_w0 > 0);  % is the bias used or not?

% This is the status of the optimization procedure.
% For GLPK, this is the exit code; see GLPKMEX for details.
% For Matlab LINPROG, if negative then no solution is found.

status = 1;


% If more than 2 classes, train the classifier one-against-all.
if C > 2,

% W1 = mclassc(a,mapping(mfilename,{is_w0,type,par,usematlab}));

  W1 = [];
  W2 = [];
  W3 = [];
  N  = [];
  for i=1:C
    mlab = 2 - (lab == i);
    aa   = dataset(+a,mlab);
    [v1,v2,v3]= fsellpc(aa,is_w0,type,par,usematlab);
    j = +v3;
    if isempty(v1),
      W1 = [];
      W2 = [];
      W3 = [];
      prwarning(1,'No solution found.');
      return;
    end
    W1 = [W1,setlabels(v1(:,1),lablist(i,:))];
    W2 = [W2;setlabels(v2(:,1),lablist(i,:))];
    W3(j) = ones(length(j),1);
    N = [N j];
  end
  [N1,N2,N3] = unique(N);
  W3 = featsel(k,N1);
  W2 = featsel(length(N1),N3)*W2;
  return

else

  Y1 = 3 - 2 * lab;   % labels     +/-1
  Y  = ones(k,1);

  alpha(1:k+1,1) = 0;

  aa = +a;
  switch upper(type)
    case {'SIMPLE'}
      f = zeros(k+z,1);
      b = -ones(m,1);
      if is_w0,
        A = -[(Y1*Y').* aa  Y1];
      else
        A = -[(Y1*Y').* aa];
      end
      [al,fval,status] = linprog(f,A,b);
      alpha(1:k+z) = al;



    case {'STANDARD'}
      L  = ones(k,1)/k;

      f  = [zeros(k+z,1); L];
      lb = [-Inf .*ones(k+z,1); zeros(k,1)];
      ub = Inf .* ones(2*k+z,1);
      b  = -ones(m,1);
      if is_w0,
        A  = -[(Y1*Y').* aa  Y1  eye(m,k)];
      else
        A  = -[(Y1*Y').* aa  eye(m,k)];
      end
      [al,fval,ststus] = linprog(f,A,b,[],[],lb,ub);
      alpha(1:k+z) = al(1:k+z);



    case {'C-SPARSE'}
      L  = ones(k,1);
      ub = Inf .* ones(3*k+z,1);
      lb = [zeros(2*k,1); -Inf.*ones(z,1); zeros(k,1)];
      b  = -ones(m,1);
      ay = (Y1*Y').* aa;
      if is_w0,
        f  = [ones(2*k,1); 0; par*L];
        A  = -[ay  -ay  Y1 eye(m,k)];
      else
        f  = [ones(2*k,1); par*L];
        A  = -[ay  -ay  eye(m,k)];
      end
      if (exist('glpkmex')>0) & (usematlab==0)
        smin     = 1;  % solve minimum
        ctype    = char(ones(m,1)*abs('U'));      % sign of inequalities
        vartype  = char(ones(3*k+z,1)*abs('C'));  % continous variables
%       lpsolver = 1;                             % Revised Simlex Method
        lpsolver = 2;                             % Interior Point Method
        params.msglev = 0; % no outputs
        [sss,hostname] = unix('hostname');
        hostname = hostname(1:end-1);
        if strcmp(hostname,'saturnus') | strcmp(hostname,'polaris') | strcmp(hostname,'neptunus')
          [al,fval,status] = glpkmex_redhat(smin,f,A,b,ctype,lb,ub,vartype,params,lpsolver);
        else
          [al,fval,status] = glpkmex(smin,f,A,b,ctype,lb,ub,vartype,params,lpsolver);
        end
      else
        [al,fval,status] = linprog (f,A,b,[],[],lb,ub);
      end
      alpha(1:k) = al(1:k) - al(k+1:2*k);
      if is_w0,
        alpha(k+1) = al(2*k+1);
      end


    case {'MU-SPARSE'}
      L   = ones(k,1)/k;
      f   = [zeros(2*k+z,1); L; -par];
      ub  = Inf .* ones(3*k+1+z,1);
      lb  = [zeros(2*k,1); -Inf.*ones(z,1); zeros(k+1,1)];
      Aeq = [ones(2*k,1); zeros(k+1+z,1)]';
      beq = 1;
      b   = zeros(m,1);
      ay  = (Y1*Y').* aa;

      if is_w0,
        A  = -[ay  -ay  Y1  eye(m,k) -ones(m,1)];
      else
        A  = -[ay -ay eye(m,k) -ones(m,1)];
      end

      if (exist('glpkmex')>0) & (usematlab==0)
        smin     = 1;  % solve minimum
        ctype    = char([ones(m,1)*abs('U'); 'S']);  % sign of inequalities
        vartype  = char(ones(3*k+1+z,1)*abs('C'));   % continous variables
%       lpsolver = 1;                                % Revised Simlex Method
        lpsolver = 2;                                % Interior Point Method
        params.msglev = 0; % no outputs
        [sss,hostname] = unix('hostname');
        hostname = hostname(1:end-1);
        if strcmp(hostname,'saturnus') | strcmp(hostname,'polaris') | strcmp(hostname,'neptunus')
          [al,fval,status] = glpkmex_redhat(smin,f,[A; Aeq],[b; beq],ctype,lb,ub,vartype,params,lpsolver);
        else
          [al,fval,status] = glpkmex(smin,f,[A; Aeq],[b; beq],ctype,lb,ub,vartype,params,lpsolver);
        end
      else
        [al,fval,status] = linprog(f,A,b,Aeq,beq,lb,ub);
      end
      alpha(1:k) = al(1:k) - al(k+1:2*k);
      if is_w0,
        alpha(k+1) = al(2*k+1);
      end

    otherwise
      disp(type)
      error ('Wrong type.');
    end
  end

  % Choose features
  ss = sum(abs(alpha(1:k)));
  J  = find(abs(alpha(1:k)) > ss*prec);

  if isempty(J) | (status <= 0) | (status > 181 | status == 150),
    prwarning(1,'Fisher classifier is trained.');
    W1 = fisherc(a);
    W2 = W1;
    W3 = featsel(k,[1:k]);
  else
    W3 = featsel(k,J);
    w = [Y; 1] .* alpha(1:k+1);
    W2 = affine(w(J),w(k+1),a(:,J),lablist,k,2);
    W2 = cnormc(W2,a(:,J));
    W1 = W3*W2;
    W1 = setname(W1,'FSELLPC');
    W2 = setname(W2,'FSELLPC');
  end
return;
