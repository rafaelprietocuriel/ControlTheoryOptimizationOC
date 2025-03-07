function b=explicitconnectiontime(mmObj)
%
% EXPLICITCONNECTIONTIME determines if the model explicitly depends on the
% connection time for multistage models.
b=[];
if isempty(mmObj)
    return
end
b=zeros(1,numberofmodels(mmObj));
for ii=1:numberofmodels(mmObj)
    if strcmp(class(mmObj.Model{ii}),'multistagemodel')
        b(ii)=explicitconnectiontime(mmObj.Model{ii});
    end
end