function idx=findresult(ocObj,reslt,varargin)
%

resultNames=fieldnames(result(ocObj));
if isempty(resultNames)
    idx=[];
    return
end

idx=strmatch(reslt,resultNames,'exact');
