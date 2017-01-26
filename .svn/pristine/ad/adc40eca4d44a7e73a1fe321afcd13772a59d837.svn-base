%SVMTRAINC Stats Support Vector Classifier (Matlab Stats Toolbox)
%
%   W = SVMTRAINC(A,KERNEL,C,OPTTYPE)
%   W = A*SVMTRAINC(KERNEL,C,OPTTYPE)
%   D = B*W
% 
% INPUT
%   A	      A PRTools dataset used fro training
%   KERNEL  Untrained mapping to compute kernel by A*(A*KERNEL) during
%           training, or B*(A*KERNEL) during testing with dataset B.
%           Default: linear kernel (PROXM('p',1))
%   C       Regularization ('boxconstraint' in SVMTRAIN)
%   OPTTYPE Desired optimizer, 'SMO' (default) or 'QP'.
%   B       PRTools dataset used for testing
% 
% OUTPUT
%   W       Mapping: Support Vector Classifier
%   D       PRTools dataset with classification results
% 
% DESCRIPTION
% This is a PRTools interface to the support vector classifier SVMTRAIN
% in Matlab's Stats toolbox. It is an alternative for STATSSVC.
%
% The evaluation of W in D = B*W makes use of the SVMCLASSIFY routine in
% the Stats toolbox. This routine outputs just labels. Consequently the
% classification matrix D has for every object the value one in the column
% corresponding with a output label and zeros for all other columns. For
% multi-class datasets the one-against-rest procedure is followed (MCLASSC)
% which may results in object classifications with zeros in all columns
% as well as multiple ones. A trained combiner may solve this, e.g.
% W = A*(SVMTRAINC*LDC)
%
% Use STATSSVC for a standard PRTools evaluation of the classifier.
%
% SEE ALSO
% DATASETS, MAPPINGS, STATSSVC, SVMTRAIN, SVMCLASSIFY, LDC

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com


function out = svmtrainc(varargin)

  checktoolbox('stats_svmtrain');
  mapname = 'StatsSVM';
  argin = shiftargin(varargin,{'prmapping','char'});
  argin = setdefaults(argin,[],proxm([],'p',1),1,'SMO');
  
  if mapping_task(argin,'definition')
    
   out = define_mapping(argin,'untrained',mapname);
   
  else
    [a,kernel,C,opttype] = deal(argin{:});
    
    if ~(ismapping(kernel) && istrained(kernel)) % training
      isdataset(a);
      islabtype(a,'crisp');
      a = testdatasize(a,'objects');
    
      % remove too small classes, escape in case no two classes are left
      [a,m,k,c,lablist,L,out] = cleandset(a,1); 
      if ~isempty(out), return; end
      
      if c > 2                        % solve multi-class case by recursion
        u = feval(mfilename,[],kernel);
        out = mclassc(a,u);           % concatenation of one-against-rest
        out = allclass(out,lablist,L);% handle with missing classes
      else                            % two class case
        labels = getlabels(a);  
        ismapping(kernel);
        isuntrained(kernel);
        prkernel(kernel);          % make kernel mapping known to prkernel
        
        pp = prrmpath('stats','svmtrain');      % check / correct path
        finishup = onCleanup(@() addpath(pp));  % restore afterwards
        if strcmpi(opttype,'SMO')
          ss = svmtrain(+a,labels,'kernel_function',@prkernel, ...
               'boxconstraint',C);
        elseif strcmpi(opttype,'QP')
          ss = svmtrain(+a,labels,'kernel_function',@prkernel, ...
               'method','qp', 'boxconstraint',C);
        else
          error('Unknown optimizer')
        end
        out = trained_mapping(a,{ss,kernel},getsize(a,3));
        %out = cnormc(out,a);       % normalise outputs for confidences
        out = setname(out,mapname);
      end
    
    else                           % evaluation
      w = kernel;                  % trained classifier
      [ss,kernel] = getdata(w);    % get datafrom training
      ismapping(kernel);
      labout = svmclassify(ss,+a); % use stats toolbox for evaluation
      nlabout = renumlab(labout,getlabels(w)); % label indices
      out = zeros(size(a,1),2);    % construct classification matrix
      out(sub2ind(size(out),[1:size(out,1)]',nlabout))= ones(size(a,1),1);
      out = setdat(a,out,w);
    end
    
  end
  
return
  



        
