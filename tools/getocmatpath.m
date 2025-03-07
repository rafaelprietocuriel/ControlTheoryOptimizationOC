function ocmatpath=getocmatpath(varargin)
% returns main directory of ocmat
%
% OCMATPATH=GETOCMATPATH() the default ocmat name is used
%
% OCMATPATH=GETOCMATPATH(OCMATNAME) the ocmat name 'OCMATNAME' is used

ocmatmainname='';
if nargin>=1
    ocmatmainname=varargin{1};
end

if isempty(ocmatmainname)
    ocmatmainname='ocmat';
end
l=length(ocmatmainname);
ocmatpath='';
pstring=path;
idx=regexp(pstring,[ocmatmainname '[;\\?]'],'once');
if isempty(idx)
    ocmatmsg(['Folder ' ocmatmainname ' is not on the MATLAB path.'])
    return
end
f=findstr(pstring(1:idx+l-1),';');

if isempty(f)
    ocmatpath=pstring(1:idx+l-1);
else
    ocmatpath=pstring(f(end)+1:idx+l-1);
end

% test if directory exists
if ~isdir(ocmatpath)
    ocmaterror('''%s'' is not a directory.',ocmatpath);
end