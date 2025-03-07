function [dUdX,dUdtsym,JUDXU]=generateexplicitcontroldynamics(m,arcid)

ctrl=retrievemodelinformation(m,'controlname');
explicitnonlinearcontrolindex=retrievemodelinformation(m,'explicitnonlinearcontrolindex',num2str(arcid));
implicitnonlinearcontrolindex=retrievemodelinformation(m,'implicitnonlinearcontrolindex',num2str(arcid));

ctrlval=control(m,[],arcid);
ai=find(constraintcombinationindex(m,arcid));
lm=lagrangemultiplier(m);
lms=lagrangemultiplier(m,[],arcid,1);
lms=lms(ai);
lmv=lm(ai);
L=lagrangian(m,[],arcid,1);
dXdt=canonicalsystem(m,[],arcid,1);
Ucell=ctrl.value(implicitnonlinearcontrolindex.value);
U=cell2vectorstring(Ucell);
% implicit control dynamics Jacobian
for ii=1:length(lmv);
    dXdt=ocmatsubs(dXdt,[lmv{ii} '=' char(lms(ii))]);
    L=ocmatsubs(L,[lmv{ii} '=' char(lms(ii))]);
end

X=cell2vectorstring([state(m) costate(m)]);
for ii=1:length(explicitnonlinearcontrolindex.value)
    idx=explicitnonlinearcontrolindex.value(ii);
    dXdt=ocmatsubs(dXdt,[ctrl.value{idx} '=' char(ctrlval(idx))]);
    L=ocmatsubs(L,[ctrl.value{idx} '=' char(ctrlval(idx))]);
end
LU=jacobian(L,U);
LUX=jacobian(LU,X);
LUU=jacobian(LU,U);
dUdtsym=-inv(LUU)*LUX*dXdt;
dUdX=-inv(LUU)*LUX;
XU=cell2vectorstring([state(m) costate(m) Ucell]);
JUDXU=(jacobian(dUdtsym,XU));