function out=exogenousstate(mmObj,solObj,varargin)
out=[];
connectflag=[];
if isempty(mmObj)
    return
end
if nargin>=3
    connectflag=varargin{1};
end
if isempty(connectflag)
    connectflag=0;
end

if ~ismmultipath(solObj)
    return
end

if numberofparts(solObj)~=numberofmodels(mmObj)
    return
end
ct=connectiontime(mmObj,solObj);

if connectflag
    for ii=1:numberofmodels(mmObj)
        out=[out exogenousstate(mmObj.Model{ii},solObj(ii),ct,1)];
    end
else
    for ii=1:numberofmodels(mmObj)
        out{ii}=exogenousstate(mmObj.Model{ii},solObj(ii),ct,connectflag,varargin{2:end});
    end
end