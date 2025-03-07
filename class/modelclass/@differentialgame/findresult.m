function idx=findresult(dgObj,reslt,varargin)
%

resultNames=fieldnames(result(dgObj));
if isempty(resultNames)
    idx=[];
    return
end

idx=strmatch(reslt,resultNames,'exact');
