function ocAsymM=heteroclinicep2ep2ocmultipath(ocObj,sol)


par=sol.modelparameter;
ocObj=changeparametervalue(ocObj,par);
arcposition=find(diff(sol.x)==0);
leftarcindex=[1 arcposition+1];
rightarcindex=[arcposition numel(sol.x)];

for ii=1:2
    hatx.y=sol.solverinfo.parameters(sol.solverinfo.equilibriumcoord{ii});
    hatx.x=0;
    hatx.arcarg=sol.solverinfo.arcarg{ii}(end);
    hatx.linearization=jacobian(ocObj,hatx);
    ocEP=dynprimitive(hatx);


    trjidx=find(sol.solverinfo.solutionindex==ii);
    trjleftarcindex=leftarcindex(trjidx);
    trjrightarcindex=rightarcindex(trjidx);

    ocTrj.x=sol.x(trjleftarcindex(1):trjrightarcindex(end))-sol.x(trjleftarcindex(1));
    ocTrj.y=sol.y(:,trjleftarcindex(1):trjrightarcindex(end));
    ocTrj.arcarg=sol.solverinfo.arcarg{ii};
    ocTrj.arcinterval=sol.solverinfo.arcinterval{ii};
    ocTrj.timehorizon=ocTrj.arcinterval(end);
    ocTrj.arcposition=[trjleftarcindex;trjrightarcindex];
    ocTrj.arcposition=ocTrj.arcposition-ocTrj.arcposition(1,1)+1;
    ocTrj.modelparameter=sol.modelparameter;
    ocTrj.modelname=sol.modelname;
    ocTrj.solverinfo.pathtype=sol.solverinfo.pathtype{ii};
    ocAsymM{ii}=ocasymptotic(octrajectory(ocTrj),ocEP);

end
ocAsymM=ocmultipath(ocAsymM);