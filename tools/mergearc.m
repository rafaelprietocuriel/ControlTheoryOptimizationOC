function sol=mergearc(sol,idx)

if nargin==1 || isempty(idx)
    idx=find(diff(sol.arcarg)==0);
end
absx=transform2absx(sol);
if ~isempty(idx)
    removex=sol.arcposition(1,idx+1);
    absx(removex)=[];
    sol.x(:,removex)=[];
    sol.y(:,removex)=[];
end
sol.arcarg(idx)=[];
sol.arcinterval(idx+1)=[];
arcposition=find(diff(sol.x)==0);
sol.arcposition=[1 arcposition+1;arcposition numel(sol.x)];
numarc=length(sol.arcposition(1,:));
for ii=1:numarc
    sol.x(sol.arcposition(1,ii):sol.arcposition(2,ii))=ii-1+transform2unit(absx(sol.arcposition(1,ii):sol.arcposition(2,ii)));
end

%--------------------------------------------------------------------------
function xnew=transform2unit(x)
xnew=[];
n=numel(x);
if n<2
    return
end
l=x(n)-x(1);
if l==0
    xnew=linspace(0,1,n);
    return
end

xnew=(x-x(1))/l;