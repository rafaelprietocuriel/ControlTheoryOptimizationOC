function out=initialstate(ocgTrj,varargin)
out=[];
if isempty(ocgTrj)
    return
end
X=state(ocgTrj,varargin{:});
if ~isempty(X)
    out=X(:,1);
end
