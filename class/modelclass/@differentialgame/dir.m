function files=dir(dgObj,placeholderidx,varargin)
%
% findmodeldata

format='';
sortflag=[];
mode='';
if isempty(dgObj)
    return
end

if nargin==1
    placeholderidx='*';
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
if isempty(placeholderidx)
    placeholderidx='*';
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
    placeholderidx=parameterindex(dgObj,placeholderidx);
end
if ~wildcardflag
    fn=filename(dgObj,format,[],placeholderidx);
else
    %fn=[modelname(dgObj) placeholderidx];
    fn=placeholderidx;
end
savedir=fullocmatfile(userdatafolder(dgObj));

files=dir(fullfile(savedir,[fn '.mat']));
switch sortflag
    case 'date'
        [dum idx]=sort([files.datenum],2,mode);
        files=files(idx);
    case 'size'
        [dum idx]=sort([files.bytes],2,mode);
        files=files(idx);
end