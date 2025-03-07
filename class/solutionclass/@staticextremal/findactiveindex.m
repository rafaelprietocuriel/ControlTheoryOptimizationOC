function out=findactiveindex(statEx,varargin)
out=[];
if isempty(statEx)
    return
end

if nargin==1
    force=1;
    opt=defaultocoptions;
end
if nargin==2
    force=varargin{1};
    opt=defaultocoptions;
elseif nargin>2
    force=varargin{1};
    opt=varargin{2};
end
if force
    tol=getocoptions(opt,'GENERAL','ZeroDeviationTolerance');
    cstr=constraint(statEx);
    lm=lagrangemultiplier(statEx);
    out=find(lm>tol | abs(lm)+abs(cstr)<tol);
else
    out=activeindex(statEx);
end