%GETILAB retrieve identifier values given identifier names or indices
%
%   LABS=GETILAB(A,REF,J)
%
% INPUT
%    A      input dataset with ilab structure (see ENABLEILAB)
%    REF    name, cell array with names or numerical indices of
%           identifier columns
%    J      indices specifying examples to be considered (optional)
%
% OUTPUT
%    VALUES identifier values. If all identifiers are numerical, matrix is
%           returned, otherwise cell array.
%
% DESCRIPTION
% GETILAB returns values of identifiers, specified by names or numerical
% indices. This function requires that the dataset supports
% named-identifiers (see ENABLEILAB). Optional vector J specifies subset
% of examples to be considered.
%
% SEE ALSO
% ENABLEILAB, REMOVEILAB, GETILAB, SETILAB, GETILABLIST

% $Id: getilab.m,v 1.4 2005/05/03 13:00:39 serguei Exp $

function labs=getilab(a,ref,J)

  isvalidilab(a);                                                    
                                                                     
  if ~exist('ref','var') | (exist('ref','var') & isempty(ref))       
    % retrieve all identifiers                                       
    ref=1:length(a.user.ilab.names);                                 
  end                                                                
                                                                     
  if ~exist('J','var')                                               
    J=1:size(a,1);                                                   
  end                                                                
                                                                     
  % get the column to be returned                                    
                                                                     
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
                                                                     
  % retrieve ident labels                                            
  labs=a.ident(J,col);                                               
  % if possible, convert the ilabs into numerical output             
  try                                                                
    labs = cell2mat(labs);                                           
    % check if the cell array contained single objects               
    if any(size(labs) ~= [length(J), length(col)]);                  
      labs = a.ident(J,col);                                         
    end                                                              
  catch                                                              
  end                                                                
                                                                     
  return                                                             
