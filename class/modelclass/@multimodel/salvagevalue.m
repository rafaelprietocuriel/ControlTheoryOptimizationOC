function out=salvagevalue(mmObj,solObj,varargin)
out=[];
if isempty(mmObj)
    return
end

if ~ismmultipath(solObj)
    return
end

if numberofparts(solObj)~=numberofmodels(mmObj)
    return
end
ct=connectiontime(mmObj,solObj);
for ii=1:numberofmodels(mmObj)
    out{ii}=salvagevalue(mmObj.Model{ii},solObj(ii),ct);
end
