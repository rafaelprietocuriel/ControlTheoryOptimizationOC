function out=matcontresult(dgObj,varargin)

out=[];
if isempty(dgObj)
    return
end
resultStruct=result(dgObj);
if isempty(resultStruct) || ~isfield(resultStruct,'MatContContinuation')
    if nargout==1
        out=[];
    end
    return
end
out=resultStruct.MatContContinuation;
if nargin==1
    return
end
if ~isnumeric(varargin{1})
    ocmaterror('Argument is not an index array.')
elseif min(varargin{1})<1 || max(varargin{1})>length(out)
    ocmaterror('Index is not in the range of the result.')
else
    out=out(varargin{1});
end
