function ocTrj=octrajectory2ocdae(ocTrj)

if isempty(ocTrj)
    return
end

if isempty(modelname(ocTrj))
    ocTrj=octrajectoy([]);
    return
end
ocObj=stdocmodel(modelname(ocTrj));
par=modelparameter(ocTrj);
ocObj=changeparametervalue(ocObj,par);
if strcmp(method(ocTrj),'dae') && isoctrajectory(ocTrj)
    ocTrj.solverinfo.statecoordinate=1:statenum(ocObj);
    ocTrj.solverinfo.costatecoordinate=statenum(ocObj)+(1:statenum(ocObj));
    ocTrj.solverinfo.controlcoordinate=ocTrj.solverinfo.costatecoordinate(end)+(1:controlnum(ocObj));
    ocTrj.solverinfo.lagrangemultipliercoordinate=ocTrj.solverinfo.controlcoordinate(end)+(1:controlconstraintnum(ocObj));
    return
end

x=state(ocObj,ocTrj,1);
l=costate(ocObj,ocTrj,1);
u=control(ocObj,ocTrj,1);
if inequalitycontrolconstraintnum(ocObj)
    lm=lagrangemultiplier(ocObj,ocTrj,1);
else
    lm=[];
end

%t=time(ocObj,ocTrj,1);
t=independentvar(ocTrj);
[t,idx]=unique(t);

X=[x;l;u;lm];
X=X(:,idx);

ocTrj=struct(ocTrj);
ocTrj.arcarg=[];
%ocTrj.arcposition=[];
ocTrj.x=t;
ocTrj.y=X;
ocTrj.solver='dae';

ocTrj.arcinterval=ocTrj.arcinterval([1 end]);

%ocTrj.solverinfo=[];
ocTrj=octrajectory(ocTrj);