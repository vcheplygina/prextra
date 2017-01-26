% PR Extra repository
% Version 0.3 6-Jan-2005
%
%Data generation
%---------------
%genriver     Generate 'bending river'-like 2D test data
%gendati      Generate random image windows
%
%Mappings
%--------
%unitm        Fixed unit mapping
%unittm       Trainable unit mapping
%
%Classifiers
%-----------
%incsvc       Incremental support vector classifier
%gatem        Gateway mapping
%emparzenc    EM Parzen Classifier for semi-supervised learning
%emc          EM algorithm for semi-supervised learning
%linprogc     Linear Program Classifier (includes sparse solutions)
%fsellpc      Low level routine for linprogc
%lessc        Least Error in Sparse Subspaces classifier (David Tax)
%adaboostc    Adaboost classifier
%rbsvc        Automatic radial basis SVC
%kannc        Fast k-nn classifier for large datasets by annquery
%dectreec     Decision tree (fast, because it's using compiled code)
%randomforestc  Random Forest (fast, because it's using compiled code)
%
%auc          error under the curve estimator
%
%Clusteranalysis
%---------------
%modeclust    Mode seeking clustering for large datasets by annquery
% 
%Image manipulation
%------------------
%impaint      Create and paint an label overlay on an image
%impatch      Create and change a polygon label overlay on an image
%neighbh      Compute the indices of neighboring pixels
%
%Feature Selection
%-----------------
%featsel2     Pairwise feature selection
%
%ilab: user-defined named identifiers
%------------------------------------
%enableilab   enable named identifiers in a dataset (ilab)
%isvalidilab  validity check for named-identifier support in a dataset
%getilablist  retrieve list of named identifiers
%addilab      add named identifier to a dataset
%getilab      retrieve identifier values given identifier names or indices
%setilab      set the values for the per-example named-identifiers in a dataset
%removeilab   remove named identifier columns
%
%
%meta: named meta data in the dataset user field
%-----------------------------------------------
%enablemeta   enable the storage of key-value pairs in a user field
%getmeta      retrieve metadata in a dataset by key names
%setmeta      add/update key-value pair(s) to a dataset
%removemeta   remove key-value pairs

