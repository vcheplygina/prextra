%ENABLEMETA  enable the storage of key-value pairs in a user field
%
%    B=ENABLEMETA(A)
%
% INPUT
%    A   input dataset
%
% OUTPUT
%    B   output dataset with initialized meta structure
%
% DESCRIPTION
% Enables the use of named key-value pairs in a dataset. These are
% stored inside the user field of a dataset.
%
% Technically, ENABLEMETA adds a field meta to the user structure.  It
% contains the cell-arrays 'keys' and 'values'. Other fields of the user
% structure than meta are not altered.
%
% SEE ALSO
% ADDMETA, REMOVEMETA, GETMETA, SETMETA

function a=enablemeta(a)
    
    u=a.user;

    meta.keys={};
    meta.values={};
    u.meta=meta;
    
    a.user=u;
    
    return    