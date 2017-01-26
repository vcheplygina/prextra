%GETILABLIST retrieve list of named identifiers
%
%   S=GETILABLIST(A)
%
% INPUT
%    A      input dataset with ilab structure (see ENABLEILAB)
%
% OUTPUT
%    S      cell-array with identifier names.
%
% DESCRIPTION
% GETILABLIST returns a list of identifiers defined in a dataset.  The
% dataset must support the user-defined named identifiers (see
% ENABLEILAB).
%
% SEE ALSO
% ENABLEILAB, REMOVEILAB, GETILAB, SETILAB, GETILABLIST

% $Id: getilablist.m,v 1.1 2005/05/02 13:05:46 pavel Exp $

function il=getilablist(a)
    
    isvalidilab(a);
    
    il=a.user.ilab.names;
    
    return;