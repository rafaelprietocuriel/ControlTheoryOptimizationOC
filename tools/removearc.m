function sol=removearc(sol,removearcidx,varargin)
%
% REMOVEARC removes an arc of an octrajectory
%
if isempty(sol) || isempty(removearcidx)
    return
end
removearcidx=sort(removearcidx);
l=length(removearcidx);
for ii=l:-1:1
    sol=removesinglearc(sol,removearcidx(ii),varargin{:});
end

function sol=removesinglearc(sol,removearcidx,varargin)

transformx=sol.x;
arcposition=find(diff(sol.x)==0);
leftarcindex=[1 arcposition+1];
rightarcindex=[arcposition numel(sol.x)];

% last x value of the previous arc
removeidx=leftarcindex(removearcidx):rightarcindex(removearcidx);
leftarcindex(removearcidx)=[];
leftarcindex(removearcidx:end)=leftarcindex(removearcidx:end)-length(removeidx);
rightarcindex(removearcidx)=[];
rightarcindex(removearcidx:end)=rightarcindex(removearcidx:end)-length(removeidx);
transformx(removeidx)=[];

% if removearcidx>1 && removearcidx<length(sol.arcinterval)
%     dx=sol.arcinterval(removearcidx+1)-sol.arcinterval(removearcidx);
% end
% sol.arcinterval(removearcidx+2:end)=sol.arcinterval(removearcidx+2:end)-dx;
sol.arcinterval(removearcidx+1)=[];
sol.arcarg(removearcidx)=[];

sol.y(:,removeidx)=[];
sol.arcposition=[leftarcindex;rightarcindex];
numarc=length(sol.arcposition(1,:));
for ii=1:numarc
    transformx(sol.arcposition(1,ii):sol.arcposition(2,ii))=ii-1+transform2unit(transformx(sol.arcposition(1,ii):sol.arcposition(2,ii)));
end
sol.x=transformx;


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