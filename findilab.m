%FINDILAB find indices of examples having a given identifier
%
%   LABS=FINDILAB(A,NAME,VALUES)
%
% INPUT
%    A        input dataset with ilab structure (see ENABLEILAB)
%    NAME     name of identifier to search
%    VALUES   value, a vector, or a cell array with identifier values 
%
% OUTPUT
%    IND      indices of examples containing identifier VALUES
%
% DESCRIPTION
% FINDILAB returns indices of examples that contain any of specified
% values.
%
% SEE ALSO
% ENABLEILAB, REMOVEILAB, findilab, SETILAB, findilabLIST

% $Id: findilab.m,v 1.1 2005/05/03 15:00:21 pavel Exp $

function ind=findilab(a,name,values)
    
    isvalidilab(a,name);
    
    t=getilab(a,name);
    
    ind=[];
    
    if isnumeric(values)
	for i=1:length(values)
	    ind=[ind; find(t==values(i))];
	end
    end
    if iscell(values)
	for i=1:length(values)
	    ind=[ind; find(t==values{i})];
	end
    end
    
    return