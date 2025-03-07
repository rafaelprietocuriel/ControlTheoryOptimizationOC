function ocObj=loadmodel(ocTrj)
%
%
mn=modelname(ocTrj);
if isempty(mn)
    ocmatmsg('Trajectory is not assigned to an ocmodel\n')
    ocObj=stdocmodel();
    return
end
ocObj=stdocmodel(mn,[],[],0);
par=modelparameter(ocTrj);
ocObj=changeparametervalue(ocObj,par);