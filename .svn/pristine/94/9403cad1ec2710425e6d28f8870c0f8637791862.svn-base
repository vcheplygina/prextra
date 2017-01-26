function out = impaint(fig,varargin)
% out = impaint(fig)
%
% Create an overlay on the current active image and allow for painting the
% overlay.
% Possible calls:
%   impaint(h,'on');               start the overlay
%   impaint(h,'off');              remove the whole overlay
%   impaint(h,classnr,[r g b]);    define a new class with colour [r g b]
% out = impaint(h,'get');          get the label image out
%
% You can leave out h, then the current active figure is used.
%
% There are also some interactive changes that you can make:
%   t    toggles impaint on and off
%  Nc    selects another class label (prefix this command with the
%         class number N)
%  Nb    sets the brush size to N (prefix this command with the brush
%         size N)
%   s    save the current labels in Matlab variable 'lab'
%
% See also  im_update


% When nothing is defined:
if nargin<2
	varargin = {'on'};
end
% Special case, because an integer can be interpreted as a handle, we
% assume that when a single input is given, it is actually a label
if nargin==1
	if isnumeric(fig)
		varargin{1} = fig;
		fig = gcf;
	end
end
if nargin<1
	fig = gcf;
end
% make sure that the first input is the figure handle:
if ~ishandle(fig)
	varargin = {fig varargin{:}};
	fig = gcf;
end

% Check if we want to get/set the active class, or (re)define a class:
if isnumeric(varargin{1})
	allUD = get(fig,'userdata');
	if ~isfield(allUD,'impaint')
		error('This image does not have an overlay defined.');
	end
	UD=allUD.impaint; % our part of UD
	if length(varargin)==1  % we want to set/get the active class
		if varargin{1}>size(UD.msk_fg,1)
			error('This class is not defined in the overlay.');
		else
			UD.currlab = varargin{1};
		end
		% and (maybe more important) get the labels to the output
		if nargout>0
			out = UD.labels;
		end
	else                    % we want to set the colours
		% We can process a whole list of colours...
		for i=1:2:(length(varargin)-1)
			[nr,clr] = checkimpaintcolor(varargin{i},varargin{i+1});
			UD.msk_fg(nr,:) = clr;
			% Definition of the active label
			UD.currlab = nr;
		end
	end
	% store everthing in the figure:
	set(fig,'UserData',[]);
	allUD.impaint=UD;
	set(fig,'UserData',allUD);
	return;
end


% Now check all the 'string' possibilities:
switch lower(varargin{1})
case 'on'

	% Set up variables:
	allUD=get(fig,'UserData'); % existing user data
	UD = [];   % Here I store the UserData

	% Put the image on screen
	UD.im_ax = findobj(get(fig,'children'),'type','axes');
	UD.im_im = findobj(get(UD.im_ax,'children'),'type','image');
	% get the size of the image
	sz = size(get(UD.im_im,'cdata'));

	%set(fig,'busyaction','cancel','DoubleBuffer','on');

	% Label storage:
	UD.labels = zeros(sz);  % background has label=0
	UD.labsize=size(UD.labels);
	
	% Mask image
	if length(sz)==2
		mskim = zeros([sz 3]);
	elseif length(sz)==3
		mskim = zeros(sz);
	else
		error('I cannot make a mask image for this image data.');
	end
	%mskim(1:50, 1:140, 1) = 1; % just to show something...
	% Definitions of the fore- and background colour
	if length(varargin)==3
	        [nr,clr] = checkimpaintcolor(varargin{2},varargin{3});
		UD.msk_fg(nr,:) = clr;
		% Definition of the active label
		UD.currlab = nr;
	else
		UD.msk_fg = [1 0 0; 0 1 0; 0 0 1];   % standard: the first class is red.
		% Definition of the active label
		UD.currlab = 1;
	end
	UD.msk_bg = [0 0 0];
	% Definition of the brush-size
	UD.brushsz = 4;
	% define character command buffer:
	UD.charbuffer='';
	UD.state='on'; % impaint is toggled on

	% make a second axis for the overlay:
	UD.pos = get(UD.im_ax,'position');
	UD.msk_ax = axes('position', UD.pos);
	set(fig,'currentaxes',UD.msk_ax);
	UD.msk_im = imagesc(mskim);
	set(UD.msk_ax,'visible','off');
	msk2 = findobj(get(UD.msk_ax,'children'),'type','image');
	set(msk2,'alphadata',0.5);

	% set the window title:
	set(fig,'Name',sprintf('imprint on, brush: %d, class: %d',UD.brushsz,UD.currlab));
	
	% to speedup coordinate transformations later:
	UD.sz = sz;
	UD.W = sz(2)./UD.pos(3);
	UD.H = sz(1)./UD.pos(4);
	UD.pos(3) = UD.pos(3)+UD.pos(1);
	UD.pos(4) = UD.pos(4)+UD.pos(2);

	% now we have to define a callback
	oldh.btndown=get(fig,'WindowButtonDownFcn');
	set(fig,'WindowButtonDownFcn','im_update(''fill'')');
	oldh.btnup=get(fig,'WindowButtonUpFcn');
	set(fig,'WindowButtonUpFcn','im_update(''stopfill'')');
	oldh.keypress=get(fig,'KeyPressFcn');
	set(fig,'KeyPressFcn','im_update(''keypress'')');	

	% preserve previously set handlers
	UD.oldhandlers=oldh;
	
	% store everthing in the figure:
	set(fig,'UserData',[]);
	allUD.impaint=UD; % set our substructure in UD
	set(fig,'UserData',allUD);
	
	out=fig;
	return
	
case 'get'
	allUD = get(fig,'userdata');
	if ~isfield(allUD,'impaint')
		error('This image does not have an overlay defined.');
	end
	out = allUD.impaint.labels;
	
	return

case 'off'
	% Remove all the traces from the image:
	allUD = get(fig,'userdata');
	if ~isfield(allUD,'impaint')
		error('This image does not have an overlay defined.');
	end
	UD=allUD.impaint;
	delete(UD.msk_ax);

	% restore old handlers
	set(fig,'WindowButtonDownFcn',UD.oldhandlers.btndown);
	set(fig,'WindowButtonUpFcn',UD.oldhandlers.btnup);
	set(fig,'KeyPressFcn',UD.oldhandlers.keypress);	
	
	rmfield(allUD,'impaint');
	set(fig,'userdata',allUD);

	return

end

return

function [nr,clr] = checkimpaintcolor(nr,clr)
if ~isnumeric(nr)
	error('The class number should be a scalar.');
end
if (length(nr)~=1)
	error('The class number should be a scalar.');
end
if ~isnumeric(clr)
	error('The color definition should be 1x3 matrix.');
end
if (length(clr)~=3)
	error('The color definition should be 1x3 matrix.');
end
if any(clr<0)
	error('Colour values should be larger than 0.');
end
if any(clr>1)
	error('Colour values should be smaller than 1.');
end
% make it a nice row vector:
clr = clr(:)';
return
