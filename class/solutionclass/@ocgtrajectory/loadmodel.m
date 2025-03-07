function ocObj=loadmodel(ocgTrj)
%
%
mn=modelname(ocgTrj);
if isempty(mn)
    ocmatmsg('Trajectory is not assigned to an ocmodel\n')
    ocObj=stdocmodel();
    return
end
ocObj=stdocmodel(mn,[],[],0);
par=modelparameter(ocgTrj);
ocObj=changeparametervalue(ocObj,par);