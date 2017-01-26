function J = neighbh(a,nbtype)
% J = neighbh(sz,nbtype)
% J = neighbh(a,nbtype)
% 
% Compute the indices of the pixels that are in the neighborhood of a
% pixel. This is done for all pixels in a whole image and returned in
% the index array J. The pixels are enumerated in this standard order
% (for an 3x5 image):
%     1   4   7  10  13
%     2   5   8  11  14
%     3   6   9  12  15
% When there is no neighbor (i.e. the pixel is on the boundary), index 0
% is returned.
% Then J becomes for a '4NN' neighborhood:
%     0  0  2  4      % pixel 1 have neighbors 2 and 4...
%     0  1  3  5      % pixel 3 has 3 neighbors..
%     0  2  0  6
%     1  0  5  7
%     2  4  6  8
%     and so on ...
%
% We can define neighborhood types nbtype:
%   '4NN'   4-nearest neighborhood
%   '8NN'   8-nearest neighborhood
%   '9NN'   9-nearest neighborhood
%  '25NN'  25-nearest neighborhood
% It is also possible to give your own definition of the neighborhood.
% You should give the indices w.r.t. the origin, for example:
%   nbtype = [-3 0;
%             +3 0;
%              0 -1];
%

if isdataset(a)
	if ~chckcoord(a)
		error('Does not look like an image dataset');
	end
	sz = getobjsize(a);
else
	if ~(size(a)==[1 2])
		error('I want to have the size of the image');
	end
	sz = a;
end
% The total number of pixels in the image:
N = prod(sz);
% The definitions of the predefined neighborhoods:
if ischar(nbtype)
	switch nbtype
	case '4NN'
		% left    up     down    right
		pos = [-1 0; 0 -1; 0  1; 1 0];
	case '8NN'
		pos = [-1 -1; -1 0; -1 1; 0 -1; 0 1; 1 -1; 1 0; 1 1];
	case '9NN'
		pos = [-1 -1; -1 0; -1 1; 0 -1; 0 0; 0 1; 1 -1; 1 0; 1 1];
	case '25NN'
		pos = [-2 -2; -2 -1; -2  0; -2  1; -2  2; -1 -2; -1 -1;...
			 -1  0; -1  1; -1  2; 0 -2; 0 -1; 0  0; 0  1; 0  2; ...
			  1 -2; 1 -1; 1  0; 1  1; 1  2; 2 -2; 2 -1; 2  0;...
			  2  1; 2  2];
	case '4star'
		%      left   up   down right
		pos = [-2 0; 0 -2; 0  2; 2 0];
	case '4cross'
		%      upleft  upright downleft downright
		pos = [-2 -2;   2 -2;   -2  2;   2 2];
	otherwise
		error('This neighborhood is not defined');
	end
else
	if size(nbtype,2)~=2
		error('The neighborhood definition should be a nx2 matrix.');
	end
	pos = nbtype;
end
% convert these neighborhood positions into changes in indices:
dx = sz(1)*pos(:,1)' + pos(:,2)';
% in horizontal limits (lower and upper):
ckhorz = -pos(:,1)';
I = find(ckhorz<0);
ckhorz(I)=ckhorz(I)+sz(2)+1;
% and in vertical limits (lower and upper):
ckvert = -pos(:,2)';
I = find(ckvert<0);
ckvert(I)=ckvert(I)+sz(1)+1;

% Setup the general parameters
D = length(dx);

% Define the pixels, and compute their coordinates:
I = (1:N)';
hcoord = floor((I-1)/sz(1))+1;
vcoord = mod(I-1,sz(1))+1;

% here the results will be stored:
J = zeros(N,D);

% and here we go:
for i=1:D
	% Here the main work of constructing the J matrix is performed,
	% it is actually just these two lines. The trick is to shift the I
	% vector up or down and store it alongside the already created matrix
	% J.

	if dx(i)>0   % look at pixels in the + direction

		J(:,i) = [I((1+dx(i)):end); zeros(dx(i),1)];

	else         % look at pixels in the - direction

		J(:,i) = [zeros(-dx(i),1); I(1:end+dx(i))];

	end
	% Unfortunately, most of the work is actually in the checking of the
	% bounds:
	% Check for the horizontal bounds:
	% look at the pixels left:
	if pos(i,1)<0
		if ckhorz(i)
			excp = find(hcoord<=ckhorz(i));
			J(excp,i) = 0;
		end
	else  % and pixels right
		if ckhorz(i)
			excp = find(hcoord>=ckhorz(i));
			J(excp,i) = 0;
		end
	end
	% Check for the vertical bounds:
	% look at the pixels *above*:
	if pos(i,2)<0
		if ckvert(i) 
			excp = find(vcoord<=ckvert(i));
			J(excp,i) = 0;
		end
	else   % otherwise the pixels below:
		if ckvert(i)
			excp = find(vcoord>=ckvert(i));
			J(excp,i) = 0;
		end
	end
end

return


function out = chckcoord(a)
% out = chckcoord(a)
%
% check if dataset a contains pixels regularly sampled from an image
% It looks a the objsize and maybe the X- and Y-coordinates.

out = 1;

sz = getobjsize(a);

% first check if we even have some 2D objsize
if length(sz)==1
	out = 0;
	return
end

featlab = getfeatlab(a);
% now check the possible X coordinate:
xcoord = strmatch('X',featlab);
if ~isempty(xcoord)
	cx = +a(:,xcoord);
	dx = diff(cx(1:sz(1)));
	if std(dx)~=0
		out = 0;
		return
	end
end
% now check the possible Y coordinate:
ycoord = strmatch('Y',featlab);
if ~isempty(ycoord)
	cy = +a(:,ycoord);
	dy = diff(cy(1:sz(1):end));
	if std(dy)~=0
		out = 0;
		return
	end
end

return
