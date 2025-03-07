function sol=changearc(sol,changearcidx,newarcid,varargin)
%
% REMOVEARC removes an arc of an octrajectory
%

if isempty(sol) || isempty(changearcidx)
    return
end
changearcidx=sort(changearcidx);
l=length(changearcidx);
for ii=l:-1:1
    sol=changesinglearc(sol,changearcidx(ii),newarcid(ii));
end

function sol=changesinglearc(sol,changearcidx,newarcid,varargin)
sol.arcarg(changearcidx)=newarcid;

sol=mergearc(sol);

function sol=mergearc(sol)

idx=find(diff(sol.arcarg)==0);
if ~isempty(idx)
    removex=sol.arcposition(1,idx+1);
    sol.x(removex)=[];
    sol.y(:,removex)=[];
end
sol.arcarg(idx)=[];
sol.arcinterval(idx+1)=[];
arcposition=find(diff(sol.x)==0);
sol.arcposition=[1 arcposition+1;arcposition numel(sol.x)];
numarc=length(sol.arcposition(1,:));
for ii=1:numarc
    sol.x(sol.arcposition(1,ii):sol.arcposition(2,ii))=ii-1+transform2unit(sol.x(sol.arcposition(1,ii):sol.arcposition(2,ii)));
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

if min(xnew)<0
    warning('Negative time values are transformed into positive values.')
    xnew=abs(xnew);
end

if any(diff(xnew)<0)
    %warning('Negative time values are transformed into positive values.')
    xnew=sort(xnew);
end