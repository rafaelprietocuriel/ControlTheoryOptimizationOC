function [ocObj,X]=matcontresult2model(ocObj,varargin)

contidx=[];
X=[];
if isempty(ocObj)
    return
end
if nargin>=2
    if isstruct(varargin{1})
        matRes=varargin{1};
    elseif isnumeric(varargin{1})
        matRes=matcontresult(ocObj,varargin{1});
    end
end
if iscell(matRes)
    matRes=matRes{1};
end
if nargin>=3
    contidx=varargin{2};
end
if isempty(matRes)
    return
end
contvar=matRes.ContinuationSolution.userinfo.varyparameterindex;
contval=matRes.ContinuationSolution.userinfo.varyparametervalue(:,contidx);
ocObj=changeparametervalue(ocObj,contvar,contval);
X=matRes.ContinuationSolution.y(:,contidx);