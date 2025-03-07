function out=contresult(dgObj,varargin)

out=[];
if isempty(dgObj)
    return
end
if nargin<=2 || (nargin==3 && isempty(varargin{2}))
    resultStruct=result(dgObj,'Continuation');
else
    resultStruct=result(dgObj,varargin{2:end});
end
if isempty(resultStruct)
    return
end
if nargin==1 || isempty(varargin{1})
    out=resultStruct;
else
    out=resultStruct(varargin{1});
end