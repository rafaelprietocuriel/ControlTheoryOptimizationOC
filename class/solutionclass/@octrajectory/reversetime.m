function ocTrjR=reversetime(ocTrj)
%
% 

ocTrjR=ocTrj;
if isempty(ocTrj)
    return
end


y=dependentvar(ocTrj);
t=independentvar(ocTrj);
arcarg=arcargument(ocTrj);
arcpos=arcposition(ocTrj);
arcintv=arcinterval(ocTrj);

ocTrjR.y=y(:,end:-1:1);
ocTrjR.x=length(arcarg)-t(end:-1:1);
ocTrjR.arcarg=arcarg(end:-1:1);
ocTrjR.arcposition=length(ocTrjR.x)-arcpos([2 1],end:-1:1)+1;
ocTrjR.timehorizon=(-1)*ocTrjR.arcinterval(end);
ocTrjR.arcinterval=arcintv(:,end:-1:1)+ocTrjR.timehorizon;
ocTrjR.solver='';

ocTrjR.solverinfo.pathtype='';