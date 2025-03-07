function out=matcontresult(ocObj,varargin)

out=[];
if isempty(ocObj)
    return
end
resultStruct=result(ocObj);
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
elseif isempty(varargin{1})
    out=[];
    return
elseif min(varargin{1})<1 || max(varargin{1})>length(out)
    ocmaterror('Index is not in the range of the result.')
else
    out=out(varargin{1});
end
