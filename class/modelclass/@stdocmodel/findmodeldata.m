function files=findmodeldata(ocObj,placeholderidx,varargin)
%
% findmodeldata

format='';
if isempty(ocObj)
    return
end

if nargin>=3
   format=varargin{1};
end
if isempty(format)
    format='%3.4g';
end
if ischar(placeholderidx)
    placeholderidx=parameterindex(ocObj,placeholderidx);
end
fn=filename(ocObj,format,[],placeholderidx);
savedir=fullocmatfile(userdatafolder(ocObj));

files=dir(fullfile(savedir,[fn '.mat']));