function sol=mergearcs(sol,idx1,idx2)

if nargin==1 || isempty(idx1) || isempty(idx2)
    idx=find(diff(sol.arcarg)==0);
    for ii=length(idx):-1:1
        sol=mergearcs(sol,idx(ii),idx(ii)+1);
    end
    return
end
absx=transform2absx(sol);
if sol.arcarg(idx1)~=sol.arcarg(idx2)
    ocmatmsg('The arcidentifiers of the arcs are not identical.')
    return
end
removex=sol.arcposition(1,idx2);
absx(removex)=[];
sol.x(:,removex)=[];
sol.y(:,removex)=[];

sol.arcarg(idx1)=[];
sol.arcinterval(idx2)=[];
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