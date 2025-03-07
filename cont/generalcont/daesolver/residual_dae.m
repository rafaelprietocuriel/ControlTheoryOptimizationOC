function res=residual_dae(t,x,dxdt,freepar,modelpar,dae)
global OCMATCONT

res=zeros(OCMATCONT.componentnumber,length(t));
for arc = 1:OCMATCONT.arcnumber
    idxmesh=1:length(t);
    rhs=dae(t,x(:,idxmesh),freepar,modelpar);
    res(OCMATCONT.zeroordercoordinate,idxmesh)=rhs(OCMATCONT.zeroordercoordinate,:);
    res(OCMATCONT.firstordercoordinate,idxmesh)=dxdt(OCMATCONT.firstordercoordinate,idxmesh)-rhs(OCMATCONT.firstordercoordinate,idxmesh);
end