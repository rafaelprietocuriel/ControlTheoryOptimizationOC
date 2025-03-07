function out=connectiontimenumber(mmObj,varargin)
out=[];
if isempty(mmObj)
    return
end

out=zeros(1,numberofmodels(mmObj));
for ii=1:numberofmodels(mmObj)
    if strcmp(class(mmObj.Model{ii}),'multistagemodel')
        out(ii)=connectiontimenumber(mmObj.Model{ii});
    end
end
