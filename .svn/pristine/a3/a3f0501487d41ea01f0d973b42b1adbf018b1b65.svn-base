%HAUSDM Hausdorff distance between datasets of image blobs
%
%	[DH,DM] = hausdm(A,B)
%
% Computes a Hausdorff distance matrix between the sets of binary images 
% A and B, or datasets containing them as features. 
% If A and B are image sets and 
% size(A) is [may,max,na]
% size(B) is [mby,mbx,nb]
% then DH is the Hausdorff distance matrix, size(DH) = [na,nb], between 
% these sets. DM is the Modified Hausdorff distance matrix.
% Preferably na <= nb (faster).
%
% See M.-P. Dubuisson and A.K. Jain, Modified Hausdorff distance for object
%     matching, Proceedings 12th IAPR International Conference on Pattern
%     Recognition (Jerusalem, October 9-13, 1994), vol. 1, IEEE, Piscataway,
%     NJ, USA,94CH3440-5, 1994, 566-568.
%
% $Id: hausdm.m,v 1.1 2005/04/25 05:54:36 duin Exp $

% Copyright: R.P.W. Duin, duin@ph.tn.tudelft.nl
% Faculty of Applied Physics, Delft University of Technology
% P.O. Box 5046, 2600 GA Delft, The Netherlands

function [dh,dm] = hausdm(A,B)

if isdataset(A) && isdataset(B)
	[dh,dm] = hausdm(data2im(A),data2im(B));
	dh = setdata(A,dh);
	dm = setdata(A,dm);
	return
end

[dummy,dummy,na] = size(A);
[dummy,dummy,nb] = size(B);
dh = zeros(na,nb);
dm = zeros(na,nb);
for i=1:na
	a = A(:,:,i); 
	J = find(any(a));
	J = [min(J):max(J)];
	K = find(any(a'));
	K = [min(K):max(K)];
	a = double(a(K,J));
	if ~isempty(a(:))
		a = bord(a,0);
	end
	ca = contourc(a,[0.5,0.5]);
	J = find(ca(1,:) == 0.5);
	ca(:,[J J+1]) =[];
	ca = ca - repmat([1.5;1.5],1,size(ca,2));
	ca = ca/max(ca(:));
	ca = ca - repmat(max(ca,[],2)/2,1,size(ca,2));
	for j = 1:nb
		b = B(:,:,j);
		J = find(any(b));
		J = [min(J):max(J)];
		K = find(any(b'));
		K = [min(K):max(K)];
		b = double(b(K,J));
		if ~isempty(b(:))
			b = bord(b,0);
		end
		cb = contourc(b,[0.5,0.5]);
		J = find(cb(1,:) == 0.5);
		cb(:,[J J+1]) =[];
		cb = cb - repmat([1.5;1.5],1,size(cb,2));
		cb = cb/max(cb(:));
		cb = cb - repmat(max(cb,[],2)/2,1,size(cb,2));
		dab = sqrt(distm(ca',cb'));
		dh(i,j) = max(max(min(dab)),max(min(dab')));
		dm(i,j) = max(mean(min(dab)),mean(min(dab')));
	end
end
	
% C = bord(A,n,m)
% Puts a border of width m (default m=1) around image A
% and gives it value n. If n = NaN: mirror image values.
function C = bord(A,n,m);
%ipcontr(0);
if nargin == 2; m=1; end
[x,y] = size(A);
if m > min(x,y)
	mm = min(x,y);
	C = bord(A,n,mm);
	C = bord(C,n,m-mm);
	return
end
if isnan(n)
   C = [A(:,m:-1:1),A,A(:,y:-1:y-m+1)];
   C = [C(m:-1:1,:);C;C(x:-1:x-m+1,:)];
else
   bx = ones(x,m)*n;
   by = ones(m,y+2*m)*n;
   C = [by;[bx,A,bx];by];
end
return