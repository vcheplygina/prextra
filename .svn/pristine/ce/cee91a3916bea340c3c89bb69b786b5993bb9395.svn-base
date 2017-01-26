%PRTOOLS_EXPORT creates and stores PRTools export
%
% PRTOOLS_EXPORT(FLAG,EXPORTDIR)
%
% Creates, if FLAG is true an export version of PRTools in
% http://prtools.tudelft.org/files/prtools.zip.
%
% If FLAG is false (default) just zip-files in EXPORTDIR are created for
% inspection and the website is not attached.
%
% EXPORTDIR should be the Matlab/prtools_export dir in the PRLab group
% directory. Default is K:\insy\PRLab\Matlab\prtools_export
%
% The intention is that this routine works on the PC as well as on the Unix
% platform. The latter, however, has not yet been tested.
%
% Save trials are prtools_export, prtools_export(0), prtools_export(0,dir)
% and inspect the prtools_export directory.
%
% This m-file should be located in PREXTRA, which should be in the path, 
% as it needs the overruled FTP commands located in PREXTRA
%
% Bob Duin 31 Aug 2016

function prtools_export(doit,exportdir)

if nargin < 2, exportdir = 'K:\insy\PRLab\Matlab\prtools_export'; end
if nargin < 1, doit = false; end

prdir = fileparts(which(mfilename));
%exportdir = fullfile(prdir,'prtools_export'); % maybe to be located elsewwhere
curdir = pwd;
[~,prdirdir] = fileparts(prdir);
if ~strcmp(prdirdir,'prextra')
  error('prtools_export should be located in the prextra directory')
end

% create temp dir
prtmp = fullfile(prdir,'tmp');
if exist(prtmp,'dir') ~= 7
  mkdir(prtmp);
end

% get rid of prtools in the path
r = which('fisherc');
if ~isempty(r)
  rmpath(r);
end

cd(prtmp);
dd = dir;
n=strmatch('prtools',{dd.name},'exact');
if ~isempty(n)
  rmdir('prtools','s');
end
    
% checkout present PRTools version
if isunix
  s = unix('svn co https://svn-mede.ewi.tudelft.nl/MM-ICT/PRTOOLS prtools');
elseif ispc
  fprintf(1,'    Tortoise will pop up, find it and click OK 2x');
  s = dos('TortoiseProc.exe /command:checkout /url:"https://svn-mede.ewi.tudelft.nl/MM-ICT/PRTOOLS" /path:prtools');
  fprintf(1,'\n');
  input('    Ready to continue?');
else 
  error('svn checkout not implemented for Mac')
end
if (s~=0)
  error('Error in running svn: could not checkout prtools')
end

% adjust version number
addpath(fullfile(pwd,'prtools'))
pver  = prtver;
fprintf(1,'    Present PRTools version is %s\n',pver{1}.Version);
pver_new = input(['    Give new version number: [' pver{1}.Version '] '],'s');
if isempty(pver_new)
  pver_new = pver{1}.Version;
end
set_version_in_Contents(pver_new);
% Commit version change
if isunix
  s = unix(['svn ci ' fullfile('prtools','Contents.m')]);
elseif ispc
  fprintf(1,'    Tortoise will pop up to commit change version number in Contents, find it and click OK 2x')
  s = dos(['TortoiseProc.exe /command:commit /path:' fullfile('prtools','Contents.m')]);
  fprintf(1,'\n');
  input('    Ready to continue?');
else 
  error(['svn commit not implemented for Mac'])
end
  
% Delete all svn files
fprintf(1,'    Removing all svn files ... \n');
if isunix
  unix('chmod -R +w prtools');
  unix('\rm -rf prtools/.svn prtools/*/.svn');
elseif ispc
  dos('rmdir /S /Q prtools\.svn');
  dos('rmdir /S /Q prtools\@prdataset\.svn');
  dos('rmdir /S /Q prtools\@prdatafile\.svn');
  dos('rmdir /S /Q prtools\@prmapping\.svn');
  dos('rmdir /S /Q prtools\private\.svn');
end

% Create export version in dir2
today = datestr(date,'ddmmmyy');
exporttemp = fullfile(prtmp,['prtex' pver_new '.' today]);
mkdir(exporttemp);
dir1 = fileparts(which('fisherc'));
dir2 = fullfile(exporttemp,'prtools');
mkdir(dir2);
copydirm(dir1,dir2); % copy main dir
copydirm(fullfile(dir1,'private'),fullfile(dir2,'private'),'all'); % copy private as it is
copydirm(fullfile(dir1,'@prdataset'),fullfile(dir2,'@prdataset'),'%'); % copy headers only
copydirm(fullfile(dir1,'@prdatafile'),fullfile(dir2,'@prdatafile'),'%'); % copy headers only
copydirm(fullfile(dir1,'@prmapping'),fullfile(dir2,'@prmapping'),'%'); % copy headers only

% Generate pcodes
fprintf(1,'    Generating pcodes ... \n');
cd(fullfile(dir2,'@prdataset')); pcode(fullfile(dir1,'@prdataset'));  % generate pcode in dir2/@dataset
cd(fullfile(dir2,'@prdatafile')); pcode(fullfile(dir1,'@prdatafile')); % generate pcode in dir2/@datafile
cd(fullfile(dir2,'@prmapping')); pcode(fullfile(dir1,'@prmapping'));  % generate pcode in dir2/@mapping

% Create zip-files
fprintf(1,'    Creating zip files ... \n');
cd(exportdir);
copyfile('License.txt',exporttemp,'f');
copyfile('Install_notes.txt',exporttemp,'f');
copyfile('Release_notes.txt',exporttemp,'f');
pfiles = ['prtex' pver_new '.' today '.zip'];
mfiles = ['prtools' pver_new '.' today '.zip'];
zip(pfiles,{'prtools','License.txt','Install_notes.txt','Release_notes.txt'},exporttemp);
zip(mfiles,[fullfile(prtmp,'prtools') '/*']);
rmdir(prtmp,'s');
  
if doit
  fprintf(1,'    FTP PRTools zip-file to website ... \n');
  % store p-files on website
  ftpsite = ftp('prtools.tudelft.nl','prtools','pvG4z1!1');
  cd(ftpsite,'httpdocs/files');
  mput(ftpsite,pfiles);
  rename(ftpsite,pfiles,'prtools.zip');
  disp('p-files have been stored on the website.')
  disp('Please check and adjust release notes and known problems on website')
  close(ftpsite)
  
  % backup p-files and m-files
  fprintf(1,'    Storing backups of mfiles and p-files ... \n');
  [s,mess] = copyfile(pfiles,fullfile(exportdir,'prtools_old_p'));
  if s
    delete(pfiles);
    disp(['p-files have been backed-up in ' exportdir])
  else
    disp(mess); 
  end
  [s,mess] = copyfile(mfiles,fullfile(exportdir,'prtools_old_m'));
  if s
    delete(mfiles);
    disp(['m-files have been backed-up in ' exportdir])
  else
    disp(mess); 
  end
end
    
cd(curdir)

	
%COPYDIRM Copy mfiles from dir to dir
function copydirm(dir1,dir2,type)
if nargin < 3
	type = 'm';
end
if exist(dir2) ~= 7
	[parentdir,subdir] = fileparts(dir2);
	mkdir(parentdir,subdir);
end
[f1,f2] = fileparts(dir2);
%disp(['Copying ' f2]);
switch type
case {'m','%'}
	if isunix
		unix(['cp ' dir1 '/*.m ' dir2]);
	else
		list = dir([dir1 '/*.m']);
	end
case {'p'}
	if isunix
		unix(['cp ' dir1 '/*.p ' dir2]);
	else
		list = dir([dir1 '/*.p']);
	end
case {'qld'}
	if isunix
		unix(['cp ' dir1 '/qld.dl* ' dir2]);
		unix(['cp ' dir1 '/qld.mex* ' dir2]);
	else
		list = dir([dir1 '/qld.dl*']);
		list = [list; dir([dir1 '/qld.mex*'])];
  end
otherwise
	if isunix
		unix(['cp ' dir1 '/* ' dir2]);
	else
		list = dir([dir1 '/*']);
	end

end
if ~isunix
	k = length(list);
	for i = 1:k
		%disp(['File copy ' list(i).name]);
		if ~list(i).isdir
			copyfile(fullfile(dir1,list(i).name),fullfile(dir2,list(i).name));
		end
	end
end

if strcmp(type,'%')
	list = dir(dir2);
	k = length(list);
	for i=1:k
		%disp(['File reduce ' list(i).name]); 
		if ~list(i).isdir
			r = gethf(fullfile(dir2,list(i).name));
			writf(fullfile(dir2,list(i).name),r);			
		end
	end
end

%LISTN List lines specified by their line number
%
% t = listn(r,n)
% Get the lines in r given by the line numbers in n.
function t = listn(r,n)
nlchar = setstr(10); % use Windows / Unix newline char
k = [0,find(r==nlchar)];
t = [];
for j = n
    t = [t,r(k(j)+1:k(j+1))];
end
return

%GETHFT Create heading file from m-file
%
%	s = gethf(file)
%
% The heading of the given m-file consisting of all starting % lines
% is isolated and returned in s.

function s = gethf(file)
[s,ns] = readf(file);
p=grep(s,'%');
if isempty(p), p = 0; end
I = find(p - [0,p(1:length(p)-1)] ~= 1);
if p(1) == 1 & isempty(I), I = length(p) + 1; end
if length(I) > 0 & I(1) > 1
	n = I(1) - 1;
	s = listn(s,1:n);
else
	s = '';
end

%GREP Get line specific lines
%
% n = grep(r,s)
% Get the numbers of all lines in the set of lines r 
% that contain s.

function k = grep(r,s)
n = [0,find(r==newline)];
m = findstr(r,s);
[i,j] = sort([n,m]);;
q = [0,j(1:length(j)-1)]-j;
k = j(find(q>0))-1;
return

function [s,n] = readf(file)
fid = fopen(file);
if fid < 0
	error(['File ' file ' not found, probably wrong directory'])
end
[s,n] = fread(fid);
if isempty(s)
  error(['File ' file ' could not be read'])
end
s = s(:)';
fclose(fid);

function writf(file,s)
filedir = fileparts(which(file));
fid = fopen(fullfile(filedir,file),'w');
if fid < 0
	error(['File ' file ' could not be opened for writing'])
end
fwrite(fid,s);
fclose(fid);

function set_version_in_Contents(prtv)
if ~isstr(prtv)
	error('Version should be given as string, e.g. ''3.2.7''')
end
s = readf('Contents.m');
n = findstr(s,newline);
r = [s(1:n(1)), '% Version ' prtv ' ' date s(n(2):end)];
writf('Contents.m',r);
disp('    Contents.m rewritten with new version number')

function myedit(file)

try 
  edit(file);
catch
  if isunix
    s = input('    Could not find editor, OK to open vi? [yes]/no: ','s');
    if isempty(s) | strcmp(s,'yes')
      unix(['vi ' file]);
    else
      disp(['    Could not open ' file ' for editing']);
    end
  else
    disp(['    Could not open ' file ' for editing']);
  end
end
