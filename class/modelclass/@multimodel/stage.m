function out=stage(mmObj)

out=[];
if isempty(mmObj)
    return
end

out=zeros(1,numberofmodels(mmObj));
for ii=1:numberofmodels(mmObj)
    if ismultistagemodel(mmObj.Model{ii})
        out(ii)=stage(mmObj.Model{ii});
    end
end
