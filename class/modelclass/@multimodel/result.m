function out=result(mmObj,varargin)
%

out=[];
if isempty(mmObj)
    return
end

if nargin==1
    out=mmObj.Result;
else
    if ~ischar(varargin{1}) || isempty(findresult(mmObj,varargin{1}))
        ocmatmsg('Input argument is not a name of the result field.\n')
    else
        out=mmObj.Result.(varargin{1});
    end
    if nargin>=3
        if ~isnumeric(varargin{2})
            ocmaterror('Argument is not an index array.')
        elseif min(varargin{2})<1 || max(varargin{2})>length(out)
            ocmaterror('Index is not in the range of the result.')
        else
            out=out(varargin{2});
        end
    end
end