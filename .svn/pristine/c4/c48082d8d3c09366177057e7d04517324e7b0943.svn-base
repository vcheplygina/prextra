%ISVALIDMETA meta validity check for datasets
%
%   I=ISVALIDMETA(A,KEYS)
%
% INPUT
%   A       input dataset
%   KEYS  (opt) name or cell array with names of required meta keys
%
% OUTPUT
%   I       logical output: the dataset is/is not valid
%
% DESCRIPTION
% This function validates if the dataset contains meta capability i.e.
% can store additional named sets of identifiers for each object.
% The metalist containing the names of ident columns is stored in the
% dataset user field. If a name or a cell array with column names is
% supplied in KEYS parameter, these are checked for existence.

function i=isvalidmeta(a,inkeys)

    if ~exist('inkeys','var'), inkeys=[]; end
    if ischar(inkeys), inkeys={inkeys}; end
    
    isdataset(a);

    i=0;
    user=a.user;
    if isfield(user,'meta')
	t=user.meta;
	if isstruct(t) & isfield(t,'keys') & iscell(t.keys) & isfield(t,'values') & iscell(t.values)
	    i=1;
	end
    end

    if (nargout == 0) & (i == 0)
	error(['Dataset with meta structure in user field expected. Use help enablemeta for details.'])
    end
    
    % validate the required keys, if supplied
    if ~isempty(inkeys)
	for j=1:length(inkeys)
	    if ~any(strcmp(user.meta.keys,inkeys{j}))
		i=0;
		if (nargout == 0) & (i == 0)
		    error(['required key ''' inkeys{j} ''' is not present in a dataset']);
		end
	    end
	end
    end
    
    
    return
	