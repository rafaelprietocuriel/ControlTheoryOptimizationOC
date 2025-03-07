function out=result(odeObj,varargin)

out=[];
if isempty(odeObj)
    return
end

if nargin==1
    out=odeObj.Result;
else
    if ~ischar(varargin{1}) || ~isfield(odeObj.Result,varargin{1})
        ocmaterror('Input argument is not a name of the result field.')
    else
        out=odeObj.Result.(varargin{1});
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