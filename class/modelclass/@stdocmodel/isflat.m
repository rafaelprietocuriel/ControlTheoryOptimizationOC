function b=isflat(ocObj,varargin)
%
% ISSPP true if limit set is of saddle type

b=zeros(1,nargin-1);
if nargin==1
    return
end
numx=statenum(ocObj);
N=parametervalue(ocObj,'N');

fx=numx/(N+1);
for ii=1:nargin-1
    x=state(ocObj,varargin{ii});
    x=reshape(x,N+1,fx);
    b(ii)=1;
    for jj=1:fx
        b(ii)=b(ii)&abs(max(diff(x(:,jj))))<=1e-8*(1+min(abs(x(:,jj))));
    end
end