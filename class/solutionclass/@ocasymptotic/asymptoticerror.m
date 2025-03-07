function out=asymptoticerror(ocAsym)

out=[];
if isempty(ocAsym)
    return
end

y=dependentvar(ocAsym.octrajectory);
ylimset=dependentvar(ocAsym.limitset);

if size(ylimset,2)>1
    return
else
    out=norm(y(:,end)-ylimset);
end
