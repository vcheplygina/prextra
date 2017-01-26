%ISVALIDILAB validity check for named-identifier support in a dataset
%
%   I=ISVALIDILAB(A,FIELDS)
%
% INPUT
%   A       input dataset
%   FIELDS  (opt) name or cell array with required identifier names 
%
% OUTPUT
%   I       logical output: the dataset is/is not valid
%
% DESCRIPTION 
% This function validates if the dataset contains named identifiers.  IN
% order to enable the use of named identifiers in a dataset, use
% ENABLEILAB.  If a name or a cell array with identifier names is supplied
% in the FIELDS parameter, these are checked for existence.
%
% SEE ALSO
% ENABLEILAB, REMOVEILAB, GETILAB, SETILAB, GETILABLIST

% $Id: isvalidilab.m,v 1.2 2005/05/02 13:05:46 pavel Exp $

function i=isvalidilab(a,fields)

    if ~exist('fields','var'), fields=[]; end
    if ischar(fields), fields={fields}; end
    
    isdataset(a);

    i=0;
    user=a.user;
    if isfield(user,'ilab')
	t=user.ilab;
	if isstruct(t) & isfield(t,'names') & iscell(t.names)
	    i=1;
	end
    end

    if (nargout == 0) & (i == 0)
	error(['Dataset with ilab structure in user field expected. Use help enableilab for details.'])
    end
    
    % validate the required fields, if supplied
    if ~isempty(fields)
	for j=1:length(fields)
	    if ~any(strcmp(user.ilab.names,fields{j}))
		i=0;
		if (nargout == 0) & (i == 0)
		    error(['required field ''' fields{j} ''' is not present in a dataset']);
		end
	    end
	end
    end
    
    
    return
	