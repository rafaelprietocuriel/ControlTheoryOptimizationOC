function files=dir(ppdeObj,placeholderidx,varargin)
%
% findmodeldata

format='';
sortflag=[];
mode='';
if isempty(ppdeObj)
    return
end

if nargin==1
    placeholderidx='';
end

if nargin>=3
   format=varargin{1};
end
if nargin>=4
   sortflag=varargin{2};
end
if nargin>=5
   mode=varargin{3};
end
if isempty(format)
    format='%3.4g';
end
if isempty(sortflag)
    sortflag='name';
end
if isempty(mode)
    mode='ascend';
end

wildcardflag=any(strfind(placeholderidx,'*'));
if ischar(placeholderidx) && ~wildcardflag
    placeholderidx=parameterindex(ppdeObj,placeholderidx);
end
if ~wildcardflag
    fn=filename(ppdeObj,format,[],placeholderidx);
else
    fn=[modelname(ppdeObj) placeholderidx];
end
savedir=fullocmatfile(userdatafolder(ppdeObj));

files=dir(fullfile(savedir,[fn '.mat']));
switch sortflag
    case 'date'
        [dum idx]=sort([files.datenum],2,mode);
        files=files(idx);
    case 'size'
        [dum idx]=sort([files.bytes],2,mode);
        files=files(idx);
end