function nrm=diffgrad(sol1,sol2,method,intmethod1d,intmethod2d)
%

nrm=[];
if nargin==3
    intmethod1d='nearest';
end
if isempty(intmethod1d)
    intmethod1d='nearest';
end
if nargin<=4
    intmethod2d=intmethod1d;
end
if isempty(intmethod2d)
    intmethod2d=intmethod1d;
end

t=unique([sol1.t,sol2.t]);
x=unique([sol1.x,sol2.x]);

sol1=devalgrad(sol1,t,x,intmethod1d,intmethod2d);
sol2=devalgrad(sol2,t,x,intmethod1d,intmethod2d);

switch lower(method)
    case 'inf'
        nrm.u=max(abs(sol1.u(:)-sol2.u(:)));
        nrm.v=max(abs(sol1.v(:)-sol2.v(:)));
    case 'l1'
        nrm.u=integrate(integrate(abs(sol1.u-sol2.u),3,t),1,x);
        nrm.v=integrate(abs(sol1.v-sol2.v),2,t);
end
nrm.method=method;

function I=integrate(Y,dim,x)

if isempty(Y)
    I=[];
    return
end
sz=size(Y);
dxsz=sz;
sz(dim)=sz(dim)-1;
dxsz(dim)=1;
dx=shiftdim(diff(x),1-dim);
dx=repmat(dx,dxsz);

switch dim
    case 1
        I=sum(reshape((Y(1:end-1,:)+Y(2:end,:))/2,sz).*dx,dim);
    case 2
        I=sum(reshape((Y(:,1:end-1,:)+Y(:,2:end,:))/2,sz).*dx,dim);
    case 3
        I=sum((Y(:,:,1:end-1)+Y(:,:,2:end))/2.*dx,dim);
end