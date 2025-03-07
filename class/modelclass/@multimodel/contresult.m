function out=contresult(mmObj,varargin)

out=[];
if isempty(mmObj)
    return
end
if nargin<=2 || (nargin==3 && isempty(varargin{2}))
    resultStruct=result(mmObj,'Continuation');
else
    resultStruct=result(mmObj,varargin{2:end});
end
if isempty(resultStruct)
    return
end
if nargin==1 || isempty(varargin{1})
    out=resultStruct;
else
    out=resultStruct(varargin{1});
end
