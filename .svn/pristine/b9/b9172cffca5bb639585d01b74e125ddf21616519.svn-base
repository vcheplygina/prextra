function lab = impatch(fig,action,varargin)
%IMPATCH Create and change a polygon label overlay on an image
%
%      OUT = IMPATCH(H,ACTION,VARARGIN)
%
% Draw some patches as overlay over an image, given by the handle H.
% Possible actions are:
%
%  Press left button to draw lines, right button to end the patch.
% >> impatch(h,'label',2,[0 0 1])    defines label 2 in blue
% >> impatch(h,'label',3,'g')        defines label 2 in blue
% >> impatch(h,'label',1)            set active label to 1
% >> out = impatch(h,'get')          get the labels out.
%  Press ESC to cancel drawing and deleting the polygon.
%
% When H is not given, it will take the current active figure.
%
% Example (when you have dipimage in your path):
%   im = laplace;
%   impaint                % now click somewhere in the image,
%                           % use leftclick to add more lines
%                           % use rightclick to close the polygon
%   impaint('label',2,'b') % make a second label, in blue
%   impaint('label',3,[0 1 0]) %make third label in green
%   out = impaint('get')   % get the label image out

% Copyright: D.M.J. Tax, D.M.J.Tax@prtools.org
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

if nargin<2
	action = [];
end
if nargin<1
	fig = gcbf;
end
% If the figure handle is not given, add it and move all other
% arguments:
if ischar(fig)
	varargin = {action varargin{:}};
	action = fig;
	fig = gcbf;
end
if ~ishandle(fig) & isnumeric(fig) % we forgot to put 'label' ???
	varargin = {fig action varargin{:}};
	action = 'label';
	fig = gcf;
end


% First check if the figure has already some patches:
if isempty(fig)
	fig = gcf;
	UD = get(fig,'userdata');
	if isempty(UD) | ~isfield(UD,'impatch')
		% OK, now we have to initialize it...
		% maybe for speedup??:
		set(fig,'busyaction','cancel','doublebuffer','on');
		UD.impatch = [];
		% initialize
		UD.impatch.mode = 'start';
		UD.impatch.curr = 0;  % current polygon
		UD.impatch.currlab = 1;  % current label
		UD.impatch.clrs = [1 0 0];
		UD.impatch.nrp = 0;   % no polygons defined
		UD.impatch.h = [];
		UD.impatch.hlab = [];
		% take care of the position and size
		% check if there is an axes object:
		hh = get(fig,'children');
		for i=length(hh)
			isax(i) = strcmp('axes',get(hh(i),'type'));
		end
		if sum(isax)==0
			error('I cannot find an axes!');
		end
		if sum(isax)>1
			error('There are several axes defined:-(');
		end
		hax = hh(find(isax));
		UD.impatch.units = get(hax,'units');
		set(hax,'units','normalized');
		UD.impatch.pos = get(hax,'position');
		him = get(hax,'children');
		% for recomputing the coordinates:
		UD.impatch.sz = size(get(him,'cdata'));
		UD.impatch.W = UD.impatch.sz(2)./UD.impatch.pos(3);
		UD.impatch.H = UD.impatch.sz(1)./UD.impatch.pos(4);
		UD.impatch.pos(3) = UD.impatch.pos(3)+UD.impatch.pos(1);
		UD.impatch.pos(4) = UD.impatch.pos(4)+UD.impatch.pos(2);

		% take care for the callback function:
		%get(fig,'windowbuttondownfcn')  %DXD
		UD.impatch.olddownfcn = get(fig,'windowbuttondownfcn');
		UD.impatch.oldmotionfcn = get(fig,'windowbuttonmotionfcn');
		UD.impatch.oldkeypressfcn = get(fig,'keypressfcn');
		set(fig,'userdata',UD);
		
		set(fig,'WindowButtonDownFcn','impatch');
		set(fig,'KeyPressFcn','impatch(''character'')');
	end
end

% The special case that we get direct actions, instead of mouse events:
if ~isempty(action)
	switch action
	case 'get'     % we want to get the label image out
		% set up:
		UD = get(gcf,'userdata');
		% go over all the patches:
		lab = zeros(UD.impatch.sz);
		[X,Y] = meshgrid(1:UD.impatch.sz(2), 1:UD.impatch.sz(1));
		for i=1:length(UD.impatch.h)
			vx = get(UD.impatch.h(i),'xdata');
			vy = get(UD.impatch.h(i),'ydata');
			im = inpolygon(X,Y,vx,vy);
			I = find(im);
			lab(I) = UD.impatch.hlab(i);
			%lab = im;
		end
		return
	case 'label'   % add a new label or change active label
		% set up:
		UD = get(gcf,'userdata');
		newnr = varargin{1};
		if newnr>size(UD.impatch.clrs,1)  % add a new label
			if length(varargin)<2
				error('I am expecting now a colour definition');
			end
			if (newnr>size(UD.impatch.clrs,1)+1)
				newnr=size(UD.impatch.clrs,1)+1;
			end
			% check if the color definition is acceptable:
			newclr = checkcolour(varargin{2});
			% and commit:
			UD.impatch.currlab = newnr;
			UD.impatch.clrs(newnr,:) = newclr;
		else                             % change the active label
			if newnr>0
				UD.impatch.currlab = newnr;
			end
		end
		% Save the results:
		set(gcf,'userdata',UD);
		return;
	case 'off'  % remove all traces from the figure
		UD = get(gcf,'userdata');
		if ~isempty(UD) & isfield(UD,'impatch')
			set(gcf,'windowbuttondownfcn',UD.impatch.olddownfcn);
			set(gcf,'windowbuttonmotionfcn',UD.impatch.oldmotionfcn);
			set(gcf,'keypressfcn',UD.impatch.oldkeypressfcn);
			delete(UD.impatch.h);
			UD = rmfield(UD,'impatch');
			if size(fieldnames(UD),1)==0, UD = []; end
			set(gcf,'userdata',UD);
			return
		end
	case 'character'   % A user pressed a key
		ch = double(get(gcf,'currentcharacter'));
		if isempty(ch), return, end
		UD = get(gcf,'userdata');
		switch ch
		case 27 % ESCAPE
%			switch UD.impatch.mode
%			case 'draw'  % cancel this complete line
				UD.impatch.hlab(UD.impatch.curr) = [];
				delete(UD.impatch.h(UD.impatch.curr));
				UD.impatch.h(UD.impatch.curr) = [];
				UD.impatch.nrp = UD.impatch.nrp-1;
				UD.impatch.curr = 0;
				UD.impatch.mode = 'start';
				set(gcf,'userdata',UD);
				set(fig,'windowbuttonmotionfcn',UD.impatch.oldmotionfcn);
%			case 'move'  % cancel the movement of the node
%			end
		end
		return

	end
end

% Get the userdata
UD = get(fig,'userdata');
% Check if the cursor position is inside the image:
set(fig,'units','normalized');
pos = get(fig,'currentpoint');
if pos(1)<UD.impatch.pos(1) | pos(1)>UD.impatch.pos(3) | ...
	pos(2)<UD.impatch.pos(2) | pos(2)>UD.impatch.pos(4)
	return
end
% Compute the position:
pos(1) = round((pos(1)-UD.impatch.pos(1))*UD.impatch.W);
pos(2) = round(UD.impatch.sz(1) - (pos(2)-UD.impatch.pos(2))*UD.impatch.H);

% Now do the actions that need the position of the cursor:
if strcmp(action,'move');
	switch UD.impatch.mode
	case 'draw'
		h = UD.impatch.h(UD.impatch.curr);
		xd = get(h,'xdata'); yd = get(h,'ydata');
		xd(end) = pos(1);    yd(end) = pos(2);
		set(h,'xdata',xd);   set(h,'ydata',yd);
	case 'move' % we want to move a node...
		minnode = UD.impatch.currnode;
		hh = UD.impatch.h(UD.impatch.curr);
		xd = get(hh,'xdata'); yd = get(hh,'ydata');
		xd(minnode) = pos(1); yd(minnode) = pos(2);
		set(hh,'xdata',xd);   set(hh,'ydata',yd);
	end
	% Save the results:
	set(fig,'userdata',UD);
	return
end

% Otherwise we have to check mouse buttons
switch UD.impatch.mode
case 'start'   % start a new line
	switch get(fig,'selectiontype')
	case 'alt'
		if length(UD.impatch.h)==0
			return;
		end
		% Find nearest node
		mind = inf; minp = 0; minnode = 0;
		for i=1:length(UD.impatch.h)
			hh = UD.impatch.h(i);
			coord = [get(hh,'xdata') get(hh,'ydata')];
			N = size(coord,1);
			diff = repmat(pos,N,1) - coord;
			[diff,nodenr] = min(sum(diff.*diff,2));
			if diff<mind
				mind = diff;
				minp = i;
				minnode = nodenr;
			end
		end
		% If we are close enough, we should move the node...
		if mind<100  % magic distance
			UD.impatch.curr = minp;
			UD.impatch.currnode = minnode;
			hh = UD.impatch.h(minp);
			xd = get(hh,'xdata'); yd = get(hh,'ydata');
			xd(minnode) = pos(1); yd(minnode) = pos(2);
			set(hh,'xdata',xd);   set(hh,'ydata',yd);
			UD.impatch.mode = 'move';
			set(fig,'WindowButtonMotionFcn','impatch(''move'')');
		end
	case 'normal'
		% Set the number of the current polygon
		UD.impatch.nrp = UD.impatch.nrp+1;
		UD.impatch.curr = UD.impatch.nrp;
		% Define the colour for this line
		thisclr = UD.impatch.clrs(UD.impatch.currlab,:);
		UD.impatch.h(UD.impatch.curr) = line([pos(1) pos(1)],[pos(2) pos(2)]);
		UD.impatch.hlab(UD.impatch.curr) = UD.impatch.currlab;
		set(UD.impatch.h(UD.impatch.curr),'color',thisclr);
		UD.impatch.mode = 'draw';
		set(fig,'WindowButtonMotionFcn','impatch(''move'')');
	end
	% Save the results:
	set(fig,'userdata',UD);
	return
case 'draw'
	switch get(fig,'selectiontype')
	case 'normal' % extend the line with a new node
		h = UD.impatch.h(UD.impatch.curr);
		xd = get(h,'xdata'); yd = get(h,'ydata');
		xd(end+1) = pos(1);    yd(end+1) = pos(2);
		set(h,'xdata',xd);   set(h,'ydata',yd);
	case 'alt'    % end the line, and fill the polygon
		h = UD.impatch.h(UD.impatch.curr);
		xd = get(h,'xdata'); yd = get(h,'ydata');
		% get rid of the last point, that was only glued to the pointer:
		xd(end) = []; yd(end) = [];
		delete(h);  % will this be OK??
		myclr = UD.impatch.clrs(UD.impatch.hlab(UD.impatch.curr),:);
		h = patch(xd,yd, myclr);
		set(h,'facealpha',0.5);
		UD.impatch.h(UD.impatch.curr) = h;
		UD.impatch.curr = 0;
		UD.impatch.mode = 'start';
		set(fig,'windowbuttonmotionfcn',UD.impatch.oldmotionfcn);
	end
case 'move'
	minnode = UD.impatch.currnode;
	hh = UD.impatch.h(UD.impatch.curr);
	xd = get(hh,'xdata'); yd = get(hh,'ydata');
	xd(minnode) = pos(1); yd(minnode) = pos(2);
	set(hh,'xdata',xd);   set(hh,'ydata',yd);
	UD.impatch.mode = 'start';
end

% Save the results:
set(fig,'userdata',UD);

return

function newclr = checkcolour(newclr)

if ischar(newclr)
	if length(newclr)>1
		error('I am expecting a single character');
	end
	defclrchar = 'rbgymc';
	defclrs = [1 0 0; 0 0 1; 0 1 0; 1 1 0; 1 0 1; 0 1 1];
	newclr = find(defclrchar==newclr);
	if isempty(newclr)
		error('Not a legal colour definition');
	end
	newclr = defclrs(newclr,:);
else
	newclr = newclr(:)';
	if length(newclr)~=3
		error('The colour definition should have R G and B');
	end
	if any(newclr<0)
		error('Colour values should be larger than 0.');
	end
	if any(newclr>1)
		error('Colour values should be smaller than 1.');
	end
end
return
