%SETILAB set the values for the per-example named-identifiers in a dataset
%
%   B=SETILAB(A,REF,VALUES,J)
%
% INPUT
%    A       input dataset with ilab structure (see MAKEILAB)
%    REF     name, cell array with names or numerical indices of
%            identifier columns
%    VALUES  identifier values to be set. If all values are numerical, 
%            matrix may be supplied, otherwise a cell array.
%    J       sample indices to be retrieved (optional)
%
% OUTPUT
%    B       output dataset
%
% DESCRIPTION
% SETILAB sets the values of identifiers, specified by names for all
% examples in a dataset (or for a subset, given by example indices J).
% VALUES may be a matrix, in case all values to be set are numerical, or
% a cell array otherwise.
%
% SEE ALSO
% ENABLEILAB, REMOVEILAB, GETILAB, SETILAB, GETILABLIST

function a=setilab(a,ref,labs,J)
    
    isvalidilab(a);

    if isempty(ref)
	% replace all identifiers
	ref=1:length(a.user.ilab.names);
    end
    
    if ~exist('J','var')
	J=1:size(a,1);
    end
    
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
 
    if isnumeric(labs)
	labs=num2cell(labs);
    end
    
    % retrieve ident labels
    temp=a.ident;
    temp(J,col)=labs;
    a.ident=temp;
    
    return
    