function b=issymmetric(ocObj,varargin)
%
% ISSPP true if limit set is of saddle type

b=zeros(1,nargin-1);
if nargin==1
    return
end
numx=statenum(ocObj);
N=parametervalue(ocObj,'N');
if numx>N+1
    return
end

for ii=1:nargin-1
    x=state(ocObj,varargin{ii});
    b(ii)=max(abs(x-x(end:-1:1)))<=1e-4;
end