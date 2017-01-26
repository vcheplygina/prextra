%REMOVEILAB remove named identifier columns
%
%   B=REMOVEILAB(A,REF)
%
% INPUT
%    A     input dataset with ilab structure (see MAKEILAB)
%    REF   name, cell array with names or numerical indices of
%          identifier columns
% OUTPUT
%    B     output dataset
%
% DESCRIPTION 
% This cunction removes named identifiers (both names and per-example
% values) from a ilab-enabled dataset (see ENABLEILAB).  If no identifier
% name or index is given (REF is empty), the named-identifier capability is
% removed from the dataset.
%
% SEE ALSO
% ENABLEILAB, REMOVEILAB, GETILAB, SETILAB, GETILABLIST

% $Id: removeilab.m,v 1.2 2005/05/02 13:05:46 pavel Exp $

function a=removeilab(a,ref)
    
    isvalidilab(a);

    if ~exist('ref','var')
	u=a.user;
	u=rmfield(u,'ilab');
	a.user=u;
	prwarning(1,'ilab capability was removed from a dataset');
	return
    end
    
    u=a.user;
    
    % get the columns to be set
    col=[];
    if isnumeric(ref)
	col=ref;    
    else
	if ischar(ref)
	    ref={ref}; % turn it into cell array
	end
	if ~iscell(ref)
	    error('numerical index, name or cell array of names expected')
	end
	% we know ref is a cell array
	for i=1:length(ref)
	    ind=strcmp(a.user.ilab.names,ref{i});
	    if sum(ind)==0
		error(['identifier label ''' ref{i} ''' not defined!']); 
	    end
	    col=[col find(ind==1)];
	end
    end
    
    if max(col)>size(a.user.ilab.names,2)
	error('index exceeds number of defined ident labels');
    end

    names=u.ilab.names;
    ncol=sort(setdiff(1:length(a.user.ilab.names),col));
    names=names(ncol);
    u.ilab.names=names;
    a.user=u;
    
    % update ident
    ident=a.ident;
    ident(:,col)=[];
    a.ident=ident;
    
    return