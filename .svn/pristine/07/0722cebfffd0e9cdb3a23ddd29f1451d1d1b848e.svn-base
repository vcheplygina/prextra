function im_update(arg1)
% im_update(arg1)
%
% Auxiliary function for impaint.m It takes care for the callbacks when
% the mouse button is clicked or released.
%
% See also  impaint.m

    fig = gcbf;
    allUD = get(fig,'userdata');
    UD=allUD.impaint; % our portion of UD

    switch lower(arg1)
     case 'fill'
      % Check if it is inside the image:
      set(fig,'units','normalized');
      pos = get(fig,'currentpoint');
      if pos(1)<UD.pos(1) | pos(1)>UD.pos(3) | pos(2)<UD.pos(2) | pos(2)>UD.pos(4)
	  return
      end
      % Compute the position:
      pos(1) = round((pos(1)-UD.pos(1))*UD.W);
      pos(2) = round(UD.sz(1) - (pos(2)-UD.pos(2))*UD.H);

      switch get(fig,'SelectionType')
       case 'normal'    % left click
			%disp('left click');
			t=floor((UD.brushsz-1)/2);
			I = -t:t;
			% just paint it foreground:
			cdata = get(UD.msk_im,'cdata');
			cdata(pos(2)+I,pos(1)+I,1) = UD.msk_fg(UD.currlab,1);
			cdata(pos(2)+I,pos(1)+I,2) = UD.msk_fg(UD.currlab,2);
			cdata(pos(2)+I,pos(1)+I,3) = UD.msk_fg(UD.currlab,3);
			set(UD.msk_im,'cdata',cdata);
			UD.labels(pos(2)+I,pos(1)+I) = UD.currlab;
       case 'alt'       % right click
			%disp('right click');
			% just paint it background:
			t=floor((UD.brushsz-1)/2);
			I = -t:t;
			cdata = get(UD.msk_im,'cdata');
			cdata(pos(2)+I,pos(1)+I,1) = UD.msk_bg(1);
			cdata(pos(2)+I,pos(1)+I,2) = UD.msk_bg(2);
			cdata(pos(2)+I,pos(1)+I,3) = UD.msk_bg(3);
			set(UD.msk_im,'cdata',cdata);
			UD.labels(pos(2)+I,pos(1)+I) = 0;
      end
      % now I want to trace the cursor:
      set(fig,'WindowButtonMotionFcn','im_update(''fill'')');
      
     case 'stopfill'
      % stop tracing the cursor:
      set(fig,'WindowButtonMotionFcn',[]);

     case 'keypress-hidden'
      t=get(fig,'currentchar');  
      if strcmp(UD.state,'off') & strcmp(t,'t')

	  UD.state='on';

	  % show the overlay axes:
	  set(UD.msk_im,'visible','on');

	  % use our impaint handlers (previous are stored)
	  set(fig,'WindowButtonDownFcn','im_update(''fill'')');
	  set(fig,'WindowButtonUpFcn','im_update(''stopfill'')');
	  set(fig,'KeyPressFcn','im_update(''keypress'')');	
	  
	  % change the window title
	  set(fig,'Name',sprintf('imprint on, brush: %d, class: %d',UD.brushsz,UD.currlab));
	  
      else
	  % invoke old handler
	  eval(UD.oldhandlers.keypress);
      end
      
     case 'keypress'
      t=get(fig,'currentchar');
      if double(t)==46
	  % first allow decimal point
	  UD.charbuffer=[UD.charbuffer '.'];
      else
	  t2=str2num(t); % is is a number?
	  if ~isempty(t2)
	      % it's a diggit -> add to the buffer
	      UD.charbuffer=[UD.charbuffer t];
	  else
	      % it's a letter command: execute
	      % get the current numerical content from buffer:
	      num=str2num(UD.charbuffer);

	      switch t
	       case 't' % toggle the impaint on and off
		if strcmp(UD.state,'on')
		    UD.state='off';

		    % hide the overlay axes:
		    set(UD.msk_im,'visible','off');
		    
		    % supply original handlers
		    set(fig,'WindowButtonDownFcn',UD.oldhandlers.btndown);
		    set(fig,'WindowButtonUpFcn',UD.oldhandlers.btnup);
		    
		    % keypress handler must allow to return back
		    set(fig,'KeyPressFcn','im_update(''keypress-hidden'')');	
		    
		    % change the window title
		    set(fig,'Name','imprint off');

		    set(fig,'UserData',[]);
		    allUD.impaint=UD;
		    set(fig,'UserData',allUD);
		    return
		    
		end
	       case 'b' % change the brush size
		if isempty(num)
		    fprintf(1,'impaint: use a numerical prefix followed by letter command!');
		end
		if num>0 
		    UD.brushsz=num;
		end
	       case 'c' % class selector
		if isempty(num)
		    fprintf(1,'impaint: use a numerical prefix followed by letter command!');
		end
		if num>0 
		    UD.currlab=num;
		end
	       case 'a' % alpha level
		if isempty(num)
		    fprintf(1,'impaint: use a numerical prefix followed by letter command!');
		end
		if num>=0 & num<=1
		    msk2 = findobj(get(UD.msk_ax,'children'),'type','image');
		    set(msk2,'alphadata',num);
		end
	       case 's' % save labels
		def='lab';
		if ~isempty(num), def=[def num2str(num)]; end
		name=inputdlg('Enter label vector name','Assign labels to a workspace variable',1,{def});
		if ~isempty(name)
		    lab=UD.labels;
		    lab=lab(1:UD.labsize(1),1:UD.labsize(2));
		    assignin('base',name{1},lab);
		end
	       otherwise
%		fprintf(1,'unknown command');
	      end
	      set(fig,'Name',sprintf('imprint on, brush: %d, class: %d',UD.brushsz,UD.currlab));	   
	      UD.charbuffer=''; % clear the buffer
	  end
      end
    end

    set(fig,'UserData',[]);
    allUD.impaint=UD;
    set(fig,'UserData',allUD);

    % DXD it does not seem to work....:
    %set(UD.msk_im,'erasemode','normal');

    return
