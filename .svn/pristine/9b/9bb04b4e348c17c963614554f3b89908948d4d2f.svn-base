%REMOVEMETA remove key-value pairs
%
%   B=REMOVEMETA(A,KEYS)
%
% INPUT
%    A     input dataset with meta structure (see ENABLEMETA)
%    KEYS  key (string) or a cell array with keys
%
% OUTPUT
%    B     output dataset
%
% DESCRIPTION
% Remove the key-value pairs specified by the keys. If no keys are
% specified, the meta facility is removed from the dataset.
%
% SEE ALSO
% REMOVEMETA, GETMETA, SETMETA

% $Id: removemeta.m,v 1.3 2005/05/02 14:12:10 pavel Exp $

function a=removemeta(a,inkeys)
    
    isvalidmeta(a);

    u=a.user;    
    if ~exist('inkeys','var')
	u=rmfield(u,'meta');
	a.user=u;
	prwarning(1,'meta capability was removed from a dataset');
	return
    end
	
    keys=u.meta.keys;
    values=u.meta.values;    
    
    if ischar(inkeys)
	inkeys={inkeys}; % turn it into cell array
    end
    if ~iscell(inkeys)
	error('name or a cell array of names expected as keys')
    end
    % we know inkeys is a cell array
    for i=1:length(inkeys)
	ind=strcmp(a.user.meta.keys,inkeys{i});
	if sum(ind)==0
	    error(['key ''' inkeys{i} ''' not defined!']); 
	end
	ind=find(ind);
	keys(ind)=[];
	values(ind)=[];
    end
    
    u.meta.keys=keys;
    u.meta.values=values;
    a.user=u;
    
    return