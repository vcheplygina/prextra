%HELPTOHTML
%
%   HELPTOHTML(TOOLBOX,DIR)
%
% Converts the help text of TOOLBOX, provided it is PRLab format, into html
% and stores them in DIR.
%
% Defaults: prtools, /data/franklin/httpd/prtools/html/prhtml

function helptohtml(toolbox,dirname)

if nargin < 2, dirname = []; end
if nargin < 1, toolbox = 'prtools'; end

actdir = cd;

switch toolbox
	case 'prtools'
		if isempty(dirname)
			if isunix
				dirname = '/data/franklin/httpd/prtools/html/prhtml';
			elseif ispc
				dirname = 's:/prt/prhtml';
			end
		end
		cd(dirname)
		prhelp2html;
		prcontent2html('prtools');
		
	otherwise
		if ~isempty(dirname)
			cd(dirname);
			prcontent2html(toolbox);
		end
end

cd(actdir)

return




