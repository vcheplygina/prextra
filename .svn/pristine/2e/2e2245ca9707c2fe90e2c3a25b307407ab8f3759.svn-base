%GETMETA retrieve metadata in a dataset by key names
%
%   [VALUES,KEYS]=GETMETA(A,KEYS)
%   [VALUES,KEYS]=GETMETA(A)
%
% INPUT
%    A       input dataset with meta structure (see MAKEMETA)
%    KEYS    key (string) or a cell array of keys (optional)
%
% OUTPUT
%    VALUES  values to be returned. Multiple values are returned in
%            a cell array.
%    KEYS    a cell-array with all keys
%
% DESCRIPTION
% Retrieve values for the given keys. If no keys are specified,
% all values are returned.
%
% SEE ALSO
% REMOVEMETA, GETMETA, SETMETA

% $Id: getmeta.m,v 1.2 2005/05/02 14:10:08 pavel Exp $

function [out,keys]=getmeta(a,inkeys)
    
    isvalidmeta(a);

    if ~exist('inkeys','var')
	out=a.user.meta.values;
	if nargout==2, keys=a.user.meta.keys; end
	return
    end
    
    u=a.user;
    keys=u.meta.keys;
    values=u.meta.values;    
    
    if ischar(inkeys)
	inkeys={inkeys}; % turn it into cell array
    end
    if ~iscell(inkeys)
	error('name or a cell array of names expected as keys')
    end
    % we know inkeys is a cell array
    out={};
    for i=1:length(inkeys)
	ind=strcmp(a.user.meta.keys,inkeys{i});
	if sum(ind)==0
	    error(['key ''' inkeys{i} ''' not defined!']); 
	end
	ind=find(ind);
	out{i}=values{ind};
    end
    
    if length(out)==1
	out=out{1};
    end

    return