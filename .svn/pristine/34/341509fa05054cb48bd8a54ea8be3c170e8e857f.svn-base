%GATEM  Gate classifier
%
%       W = GATEM(A,WGATE,W)
%
% Make a gate classifier WGATE that splits the feature space in subsets,
% train a classifier W{i} on each of the subsets and combine the results
% again. Thus W is a cell-array (cell-vector) of untrained classifiers.
% Both the gate classifier as the secondary classifiers are trained on
% dataset A.  This gatem is like a mixture of experts, where each
% classifier is focussing on a special area in the feature space.
%
% Example:  a = gendatb;    % generate a 2-class banana dataset
%           w = gatem(a, ldc, {parzenc, qdc});
%                           % split the feature space in two by an ldc
%                           % train a parzenc in one half, and a qdc in
%                           % the other half
%           scatterd(a); plotc(w)
%
% In this version the ordering of the data that is passed to the
% classifiers W{i} is based on the ordering as defined by
% getlabels(A*WGATE).  Probably this should be made more flexible in a
% future version.

% Copyright: D.M.J. Tax, D.M.J.Tax@prtools.org
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function w = gatem(a,wgate,W)

if nargin < 1 | isempty(a) 
	% When no inputs are given, we are expected to return an empty
	% mapping:
	W = mapping(mfilename,{wgate,W});
	W = setname(W,'Gate classifier');
	return
end

if ~strcmp(getmapping_file(wgate),mfilename)  % training
	% training consist of training the gate classifier, splitting the
	% dataset A according to the output label, and training the
	% individual classifiers on the smaller datasets.

	% First train the gate:
	if ~istrained(wgate)
		wgate = a*wgate;
	end
	gatelabels = getlabels(wgate);
	[k,c] = size(wgate);

	% Do some checking
	if length(W)~=c
		error('The number of classifiers in W should be equal to the number of classes');
	end

	% Now map the data through the classifier:
	out = a*wgate;
	% and extract the indices of the sub-datasets:
	[dummy,I] = max(+out,[],2);

	% Train the secondary mappings:
	for i=1:c
		%Empty classifiers and trained classifiers don't have to be trained
		if ~isempty(W{i}) & ~istrained(W{i}) 
			% extract the datasets,
			b = a(find(I==i),:);
			% be careful, this can go wrong:
			if isempty(b)
				warning('One of the classifiers did not get any training data');
			else
				% and train the mapping
				W{i} = b*W{i};
			end
		end
		
		% Does this classifier output the same as the gate??
		if ~isempty(W{i})
			labels{i} = matchlablist(getlabels(W{i}),gatelabels);
		else
			labels{i} = [];
		end
	end

	% Now store everything:
	V.wgate = wgate;
	V.W = W;
	V.labels = labels;
	w = mapping(mfilename,'trained',V,gatelabels,k,c);
	w = setname(w,'Gate classifier');

else                  % evaluation
	% get the data out:
	V = getdata(wgate);
	[n,k] = size(a);

	% Apply the gate classifier:
	out = +(a*V.wgate);
	[dummy,lab] = max(+out,[],2);

	% Now apply the secondary classifiers
	for i=1:length(V.W)
		if ~isempty(V.W{i})
			% extract the datasets,
			I = find(lab==i);
			if ~isempty(I)  % so, if there is data to classify:
				% compute the output:
				tmpout = +(a(I,:)*V.W{i});
				% store the output:
				%out(I,:) = tmpout(:,V.labels{i});
				out(I,V.labels{i}) = tmpout;
			end
		end
	end

	% Fill in the data, keeping all other fields in the dataset intact
	w = setdata(a,out);
	w = set(w,'featlab',getlabels(wgate),'featsize',getsize_out(wgate));

end

return
