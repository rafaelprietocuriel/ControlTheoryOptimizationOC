function sol=generateodestruct(dynPrim,initnummesh,timesettransformation)

sol=[];
if isempty(dynPrim)
    return
end

if isequilibrium(dynPrim)
    if timesettransformation.asymptoticapproximation>=1
        sol.x=[0 logspace(log10(1/timesettransformation.asymptoticapproximation),0,initnummesh-1)];
    else
        sol.x=linspace(0,1,initnummesh);
    end
    sol.y=repmat(dependentvar(dynPrim),1,initnummesh);
    sol.parameters=[];
    sol.arcinterval=[0 timesettransformation.asymptoticapproximation]; %
    sol.timehorizon=inf;
end
sol.arcarg=arcargument(dynPrim);
sol.x0=sol.x(1);
sol.idata.tangent=[];
sol.idata.coeff=[];