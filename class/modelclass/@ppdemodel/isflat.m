function b=isflat(ppdeObj,varargin)
%

b=zeros(1,nargin-1);
if nargin==1
    return
end
for ii=1:nargin-1
    x=state(ppdeObj,varargin{ii});
    b(ii)=abs(max(diff(x)))<=1e-8*min(abs(x));
end