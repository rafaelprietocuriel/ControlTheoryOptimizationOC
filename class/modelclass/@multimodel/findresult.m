function idx=findresult(mmObj,reslt,varargin)
%

resultNames=fieldnames(result(mmObj));
if isempty(resultNames)
    idx=[];
    return
end

idx=strmatch(reslt,resultNames,'exact');
