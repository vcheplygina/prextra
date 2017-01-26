%ADDILAB add named identifier to a dataset
%
%   B=ADDILAB(A,REF,LABS)
%
% INPUT
%    A     input dataset with ilab structure (see MAKEILAB)
%    REF   name, cell array with names or numerical indices of
%          identifier columns
%    LABS  identifiers to be set. If all identifiers are numerical, 
%          matrix may be supplied, otherwise a cell array.
%
% OUTPUT
%    B     output dataset
%
% DESCRIPTION
% This function adds one or more named identifiers to the dateset
% and sets their per-example values. To be used, the named-identifier 
% facility must be enabled in a dataset (see ENABLEILAB).
%
% SEE ALSO
% ENABLEILAB, REMOVEILAB, GETILAB, SETILAB, GETILABLIST

% $Id: addilab.m,v 1.2 2005/05/02 13:05:46 pavel Exp $

function a=addilab(a,name,labs)
    
    isvalidilab(a);

    u=a.user;
    
    % is the name already defined?
    ind=strcmp(u.ilab.names,name);
    if any(ind)
	error([name ' already defined in ilab!']);
    else
	
	% check the validity of labs to be added
	if length(labs)~=size(a,1)
	    length(labs)
	    error('wrong size of ident labels to be included!');
	end

%	labs=reshape(labs,size(a,1),1);
	
	if isnumeric(labs)
	    labs=num2cell(labs);
	end
	
	if ~iscell(labs)
	    error('ident labels to include must be either cell array or numerical vector');
	end
	
	% update ilab
	u.ilab.names{end+1}=name;
	a.user=u;
	ident=a.ident;
	ident(:,end+1)=labs;
	a.ident=ident;
	
    end
    
    return