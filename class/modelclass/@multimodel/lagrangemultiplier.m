function out=lagrangemultiplier(mmObj,solObj,varargin)
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
if connectflag
    for ii=1:numberofmodels(mmObj)
        out=[out lagrangemultiplier(mmObj.Model{ii},solObj(ii),1)];
    end
else
    out=lagrangemultiplier(mmObj.Model{1},solObj(1),1);
end