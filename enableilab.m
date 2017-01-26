%ENABLEILAB - enable named identifiers in a dataset (ilab)
%
%    B=ENABLEILAB(A)
%
% INPUT
%    A   input dataset
%
% OUTPUT
%    B   output dataset with initialized ilab structure
%
% DESCRIPTION
% Enables the use of multiple named identifiers in a dataset. For each
% example in a dataset, multiple identifiers may be defined. Each
% identifier has a name. The names (identifier labels i.e. ilab) are stored
% in a ilablist in a user field of a dataset. This function removes all the
% identifiers stored in a dataset, add one column, called ident, and
% populates it with unique indices.
%
% Technically, ENABLEILAB adds a field ilab to the user structure.  The
% field names of the ilab stores the column names of the dataset ident
% cell-array. Other fields of user structure than ilab are not altered.
%
% SEE ALSO
% ADDILAB, REMOVEILAB, GETILAB, SETILAB, GETILABLIST

function a=enableilab(a)

    % reseting the identifiers
    a.ident=num2cell((1:size(a,1))');
    prwarning(3,'identifiers were reset in the dataset');
    
    ilab.names={'ident'};
    u.ilab=ilab;
    
    a.user=u;
    
    return
