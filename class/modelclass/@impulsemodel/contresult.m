function out=contresult(ocObj,varargin)

out=[];
if isempty(ocObj)
    return
end
resultStruct=result(ocObj);
if isempty(resultStruct) || ~isfield(resultStruct,'Continuation')
    if nargout==1
        out=[];
    end
    return
end
if nargin==1
    out=resultStruct.Continuation;
else
    out=resultStruct.Continuation(varargin{1});
end
