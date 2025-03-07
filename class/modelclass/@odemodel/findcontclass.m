function idx=findcontclass(odeObj,contclass,varargin)
%

contresultStruct=contresult(odeObj);
if isempty(contresultStruct)
    idx=[];
    return
end

idx=zeros(1,length(contresultStruct));
for ii=1:length(contresultStruct)
    idx(ii)=strcmp(contclass,contresultStruct{ii}.ContinuationClassification);
end
idx=find(idx);