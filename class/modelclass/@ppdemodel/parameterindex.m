function idx=parameterindex(ppdeocObj,parvar)
idx=[];

if isempty(ppdeocObj)
    return
end
info=retrievemodelinformation(ppdeocObj.Model,'parametername');
if ischar(parvar)
    parvar=regexp(parvar, ',','split');
else
    parvar=cellstr(info.value(parvar,:));
end

numparvar=numel(parvar);
idx=zeros(1,numparvar);
for ii=1:numparvar
    ff=find(strcmp(info.value,parvar{ii}));
    if ~isempty(ff)
        idx(ii)=ff;
    end
end
idx(~idx)=[];