%SETMETA add/update key-value pair(s) to a dataset
%
%   B=SETMETA(A,KEYS,VALUES)
%
% INPUT
%    A       input dataset with meta structure (see MAKEMETA)
%    KEYS    string key or a cell-array with key names
%    VALUES  key-specific values. If more keys are defined, VALUES
%            must be stored in a cell-array
%
% OUTPUT
%    B     output dataset
%
% DESCRIPTION
% This function sets the key-value pair(s) in a user field of a dataset.
%
% SEE ALSO
% REMOVEMETA, GETMETA, SETMETA

function a=setmeta(a,inkeys,invalues)
    
    isvalidmeta(a);

    u=a.user;

    if ischar(inkeys), inkeys={inkeys}; end
    if ~iscell(invalues), invalues={invalues}; end    

    if ~iscell(inkeys) | ~iscell(invalues)
	error('unsupported format of keys or values')
    end
    if length(inkeys)~=length(invalues)
	error('keys and values to be set must have equal sizes');
    end

    % update meta in user
    u=a.user;
    for j=1:length(inkeys)
	% is the key already defined?
	ind=strcmp(u.meta.keys,inkeys{j});
	if any(ind)
	    k=find(ind);
	    u.meta.keys{k}=inkeys{j};
	    u.meta.values{k}=invalues{j};
	else % not found, add a new key-value
	    u.meta.keys{end+1}=inkeys{j};
	    u.meta.values{end+1}=invalues{j};
	end	
    end
    a.user=u;
    
    return