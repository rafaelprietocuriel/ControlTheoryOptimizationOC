function ocAsymM=heteroclinicep2lp2ocmultipath(ocObj,sol)


par=sol.modelparameter;
ocObj=changeparametervalue(ocObj,par);
arcposition=find(diff(sol.x)==0);
leftarcindex=[1 arcposition+1];
rightarcindex=[arcposition numel(sol.x)];


hatx.y=sol.solverinfo.parameters(sol.solverinfo.equilibriumcoordinate{2});
hatx.x=0;
hatx.arcarg=sol.solverinfo.arcarg{3}(end);
hatx.linearization=jacobian(ocObj,hatx);
ocEP=dynprimitive(hatx);


trjidx=find(sol.solverinfo.solutionindex==1);
trjleftarcindex=leftarcindex(trjidx);
trjrightarcindex=rightarcindex(trjidx);

ocTrj.x=sol.x(trjleftarcindex(1):trjrightarcindex(end))-sol.x(trjleftarcindex(1));
ocTrj.y=sol.y(:,trjleftarcindex(1):trjrightarcindex(end));
ocTrj.arcarg=sol.solverinfo.arcarg{1};
ocTrj.arcinterval=sol.solverinfo.arcinterval{1};
ocTrj.timehorizon=sol.solverinfo.timehorizon(1);
ocTrj.arcposition=[trjleftarcindex;trjrightarcindex];
ocTrj.arcposition=ocTrj.arcposition-ocTrj.arcposition(1,1)+1;
ocTrj.modelparameter=sol.modelparameter;
ocTrj.modelname=sol.modelname;
ocTrj.solverinfo.pathtype=sol.solverinfo.pathtype{1};

lcidx=find(sol.solverinfo.limitcycleindex==1);
lcleftarcindex=leftarcindex(lcidx);
lcrightarcindex=rightarcindex(lcidx);

ocLC.octrajectory.x=sol.x(lcleftarcindex(1):lcrightarcindex(end))-sol.x(lcleftarcindex(1));
ocLC.octrajectory.y=sol.y(:,lcleftarcindex(1):lcrightarcindex(end));
ocLC.octrajectory.arcarg=sol.solverinfo.arcarg{2};
ocLC.octrajectory.arcinterval=sol.solverinfo.arcinterval{2};
ocLC.octrajectory.timehorizon=sol.solverinfo.timehorizon(2);
ocLC.octrajectory.arcposition=[lcleftarcindex;lcrightarcindex];
ocLC.octrajectory.arcposition=ocLC.octrajectory.arcposition-ocLC.octrajectory.arcposition(1,1)+1;
ocLC.octrajectory.modelparameter=sol.modelparameter;
ocLC.octrajectory.modelname=sol.modelname;
ocLC.octrajectory.linearization=sol.solverinfo.monodromy;
ocLC.period=sol.solverinfo.timehorizon(2);

ocAsym1=ocasymptotic(octrajectory(ocTrj),ocLC);


trjidx=find(sol.solverinfo.solutionindex==3);
trjleftarcindex=leftarcindex(trjidx);
trjrightarcindex=rightarcindex(trjidx);

ocTrj.x=sol.x(trjleftarcindex(1):trjrightarcindex(end))-sol.x(trjleftarcindex(1));
ocTrj.y=sol.y(:,trjleftarcindex(1):trjrightarcindex(end));
ocTrj.arcarg=sol.solverinfo.arcarg{3};
ocTrj.arcinterval=sol.solverinfo.arcinterval{3};
ocTrj.timehorizon=sol.solverinfo.timehorizon(3);
ocTrj.arcposition=[trjleftarcindex;trjrightarcindex];
ocTrj.arcposition=ocTrj.arcposition-ocTrj.arcposition(1,1)+1;
ocTrj.modelparameter=sol.modelparameter;
ocTrj.modelname=sol.modelname;
ocTrj.solverinfo.pathtype=sol.solverinfo.pathtype{2};

ocAsym2=ocasymptotic(octrajectory(ocTrj),ocEP);
ocAsymM{1}=ocAsym1;
ocAsymM{2}=ocAsym2;
ocAsymM=ocmultipath(ocAsymM);