function [dUdX,dUdtsym,JUDX,JUDU]=generateimplicitcontroldynamics(m,arcid)

ctrl=retrievemodelinformation(m,'controlname');
explicitnonlinearcontrolindex=retrievemodelinformation(m,'explicitnonlinearcontrolindex',num2str(arcid));
implicitnonlinearcontrolindex=retrievemodelinformation(m,'implicitnonlinearcontrolindex',num2str(arcid));
ctrlval=control(m,[],arcid);
ai=find(constraintcombinationindex(m,arcid));
lm=lagrangemultiplier(m);
lms=lagrangemultiplier(m,[],arcid,1);
lms=lms(ai);
lmv=lm(ai);
Ucell=ctrl.value(implicitnonlinearcontrolindex.value);
dUdt=sym(zeros(length(Ucell),1));
for ii=1:length(Ucell)
    dUdt(ii)=sym(['D' Ucell{ii} 'Dt']);
end
U=cell2vectorstring(Ucell);

X=cell2vectorstring([state(m) costate(m)]);
L=lagrangian(m,[],arcid);
LU=jacobian(L,U);
dXdt=canonicalsystem(m,[],arcid,1);
for ii=1:length(explicitnonlinearcontrolindex.value)
    idx=explicitnonlinearcontrolindex.value(ii);
    dXdt=ocmatsubs(dXdt,[ctrl.value{idx} '=' char(ctrlval(idx))]);
    LU=ocmatsubs(LU,[ctrl.value{idx} '=' char(ctrlval(idx))]);
end

% implicit control dynamics Jacobian
for ii=1:length(lmv);
    dXdt=ocmatsubs(dXdt,[lmv{ii} '=' char(lms(ii))]);
    LU=ocmatsubs(LU,[lmv{ii} '=' char(lms(ii))]);
end
LUX=jacobian(LU,X);
LUU=jacobian(LU,U);
dUdtsym=-inv(LUU)*LUX*dXdt;
G=LUX*dXdt+LUU*dUdt;
GX=jacobian(G,X);
GU=jacobian(G,U);
dUdX=-inv(LUU)*LUX;
JUD=-inv(LUU)*(GX+GU*dUdX);
JUDU=-inv(LUU)*GU;
JUDX=simple(JUD-JUDU*dUdX);




for ii=1:length(Ucell);
    JUDX=ocmatsubs(JUDX,[char(dUdt(ii)) '=' char(dUdtsym(ii))]);
    JUDU=ocmatsubs(JUDU,[char(dUdt(ii)) '=' char(dUdtsym(ii))]);
end