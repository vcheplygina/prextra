function impaintbrush(h,sz)
% impaintbrush(h,sz)
%
% Set the brush size as it is used in impaint.m.

allUD = get(h,'UserData');

if ~isfield(allUD,'impaint')
	error('No overlay is defined.');
end
UD = allUD.impaint;

UD.brushsz = sz;

allUD.impaint = UD;
set(h,'UserData',allUD);

return


