function sol=generateodestruct(ppdePrim,initnummesh,truncationtime)

sol=[];
if isempty(ppdePrim)
    return
end

%sol.x=linspace(0,truncationtime,initnummesh);
sol.x=linspace(0,1,initnummesh);
sol.y=repmat(dependentvar(ppdePrim),1,initnummesh);
sol.parameters=[];
%sol.arcinterval=[0 truncationtime]; %
sol.arcinterval=[0 truncationtime]; %
sol.timehorizon=inf;
sol.arcarg=0;%arcargument(ppdePrim);
sol.x0=sol.x(1);
sol.idata.tangent=[];
sol.idata.coeff=[];