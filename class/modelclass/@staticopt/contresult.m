function out=contresult(ocObj,varargin)

out=[];
if isempty(ocObj)
    return
end
if nargin<=2 || (nargin==3 && isempty(varargin{2}))
    resultStruct=result(ocObj,'Continuation');
else
    resultStruct=result(ocObj,varargin{2:end});
end
if isempty(resultStruct)
    return
end
if nargin==1 || isempty(varargin{1})
    out=resultStruct;
else
    out=resultStruct(varargin{1});
end
