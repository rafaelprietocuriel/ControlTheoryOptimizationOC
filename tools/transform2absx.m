function absx=transform2absx(sol)
% transform x (relative/normalized time) into absolute time
arcposition=find(diff(sol.x)==0);
leftarcindex=[1 arcposition+1];
rightarcindex=[arcposition numel(sol.x)];

absx=sol.x;
for ii=1:numel(leftarcindex)
    xint=sol.x(leftarcindex(ii):rightarcindex(ii));
    absx(leftarcindex(ii):rightarcindex(ii))=(xint-xint(1))*(sol.arcinterval(ii+1)-sol.arcinterval(ii))+sol.arcinterval(ii);
end