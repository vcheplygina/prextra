%ECOC Error-correcting output code
% 
%      W = ECOC(A,BASECL,CM,RULE)
%
% INPUT
%    A       Dataset
%    BASECL  Untrained base classifier
%    CM      Coding matrix
%    RULE    Combination rule
%
% OUTPUT
%    W       ECOC  classifier
%
% DESCRIPTION
% Computation of the classifier using Error-correcting output codes to
% create a multi-class classifier from a base classifier/dichotomizer.
% The classes in dataset A are relabeled according to the coding matrix
% CM (containing +1, -1 or 0). The default coding matrix is a 1-vs-rest
% coding. For a three-class problem this becomes:
%     CM = [ 1 -1 -1;
%           -1  1 -1;
%           -1 -1  1];
%
% For the evaluation of objects, the outputs of the dichotomizers are
% combined using the combination rule RULE. Currently the following
% rules are implemented:
%     'none'     take the maximum output (this only works for 1-vs-all
%                coding matrices)
%     'hamming'  standard Hamming distance between the discretized
%                classifier outcomes and the coding matrix. When a zero
%                entry in the coding matrix appears, this entry is not
%                used in the distance computation
%
% NOTE: the order of the classes as they are used in the coding matrix
% is determined by the lablist of A.
%
% SEE ALSO
%    MCLASSC, MAXC

% Copyright: D.M.J. Tax, D.M.J.Tax@prtools.org
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands
  
function W = ecoc(a,basecl,CM,rule)

if nargin < 4 | isempty(rule), rule = 'hamming'; end
if nargin < 3 | isempty(CM), CM = []; end
if nargin < 2 | isempty(basecl), basecl = fisherc; end
if nargin < 1 | isempty(a) 
	W = mapping(mfilename,{w,L,rule});
	W = setname(W,'ECOC');
	return
end

if ~ismapping(basecl) | ~istrained(basecl)   %training

	[n,k,c] = getsize(a);
	nlab = getnlab(a);
	% define the coding matrix if not given:
	if isempty(CM)
		CM = repmat(-1,c,c) + 2*eye(c);
	else
		% check if the coding matrix is OK
		if size(CM,1)~=c
			error('Coding matrix does not have correct number of rows.');
		end
		if (max(CM(:))>1)|(min(CM(:))<-1)
			error('Entries in the coding matrix should be -1 or +1.');
		end
	end
	l = size(CM,2);
   w = cell(1,l);

	% Train it:
	for i=1:l
		% create a new dataset
		lab = nlab;
		for j=1:c
			lab(lab==j) = CM(j,i);
		end
		% should I remove the objects with 0 label?
		b = setlabels(a,lab);
		% train the base dichotomizer:
		w{i} = b*basecl;
	end

	%and save all useful data:
	W.CM = CM;
	W.w = w;
	W.rule = rule;
	W = mapping(mfilename,'trained',W,getlablist(a),k,c);
	W = setname(W,'ECOC');

else                               %testing

	% Extract the data:
	W = getdata(basecl);
	m = size(a,1);
	[c,l] = size(W.CM);

	% Apply all the base classifiers:
	z = zeros(m,l);
	for i=1:l
		tmp = a*W.w{i};
		id = findfeatlab(tmp,1);
		z(:,i) = +tmp(:,id);
	end
	z = 2*z-1;  %DXD classifier output should be between -1 and +1!!!
%[W.CM; z]
	% and apply the combination rule:
	out = zeros(m,c);
	switch W.rule
	case 'none'
		if size(z,2)~=c
			error('You have to combine when you don''t do one-vs-rest.');
		end
		out = z;
	case 'hamming'
		for i=1:m
			zz = sign(z(i,:)); % discrete output
			% hamming distance of object i to all code words:
			% remove entries for which the coding matrix has a 0:
			CM0 = (W.CM~=0);
			d = sum(CM0.*abs(repmat(zz,c,1)-W.CM),2);
%keyboard
			[md,Icl] = min(d);
			% class Icl won!
			out(i,Icl) = 1;
		end
	otherwise
		error('This combination rule is not implemented yet.');
	end

	% Store the outputs
	W = setdat(a,out,basecl);
end
return


