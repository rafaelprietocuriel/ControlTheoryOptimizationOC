function idx=parameterindex(odeObj,parvar)
idx=[];

if isempty(odeObj)
    return
end
info=retrievemodelinformation(odeObj.Model,'parametername');
if ischar(parvar)
    parvar=regexp(parvar, ',','split');
elseif isnumeric(parvar)
    parvar=cellstr(info.value(parvar,:));
elseif isa(parvar,'sym')
    idx=parameterindex(ocObj,char(parvar));
    return
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