function sol=multiplyentry(sol,multiplyposition,multiple)


absx=transform2absx(sol);
totalpos=sort([1:length(sol.x) multiplyposition(ones(1,multiple))]);
absx=absx(totalpos);
sol.x=sol.x(totalpos);
sol.y=sol.y(:,totalpos);
arcposition=find(diff(sol.x)==0);
leftarcindex=[1 arcposition+1];
rightarcindex=[arcposition numel(sol.x)];
numarc=numel(rightarcindex);
sol.arcinterval=[absx(leftarcindex) absx(rightarcindex(end))];
sol.arcposition=[leftarcindex;rightarcindex];
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