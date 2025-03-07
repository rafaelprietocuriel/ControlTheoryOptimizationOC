function ocAsymM=heteroclinicep2lc2ocmultipath(ocObj,sol)


par=sol.modelparameter;
ocObj=changeparametervalue(ocObj,par);
arcposition=find(diff(sol.x)==0);
leftarcindex=[1 arcposition+1];
rightarcindex=[arcposition numel(sol.x)];

hetorder=length(sol.solverinfo.pathtype);
ocAsymM=cell(1,hetorder);
for ii=1:hetorder
    trjidx=find(sol.solverinfo.solutionindex==ii);
    trjleftarcindex=leftarcindex(trjidx);
    trjrightarcindex=rightarcindex(trjidx);
    ocTrj.x=sol.x(trjleftarcindex(1):trjrightarcindex(end))-sol.x(trjleftarcindex(1));
    ocTrj.y=sol.y(:,trjleftarcindex(1):trjrightarcindex(end));
    ocTrj.arcarg=sol.solverinfo.arcarg{ii};
    ocTrj.arcinterval=sol.solverinfo.arcinterval{ii};
    ocTrj.timehorizon=sol.solverinfo.timehorizon(ii);
    ocTrj.arcposition=[trjleftarcindex;trjrightarcindex];
    ocTrj.arcposition=ocTrj.arcposition-ocTrj.arcposition(1,1)+1;
    ocTrj.modelparameter=sol.modelparameter;
    ocTrj.modelname=sol.modelname;
    ocTrj.solverinfo.pathtype=sol.solverinfo.pathtype{ii};
    switch sol.solverinfo.limitsettype{sol.solverinfo.octrajectory2limset(ii,2)}
        case 'e'
            hatx.y=sol.solverinfo.parameters(sol.solverinfo.equilibriumcoordinate{sol.solverinfo.octrajectory2limset(ii,3)});
            hatx.x=0;
            hatx.arcarg=sol.solverinfo.arcarg{ii}(end);
            hatx.linearization=jacobian(ocObj,hatx);
            hatx.modelparameter=sol.modelparameter;
            hatx.modelname=sol.modelname;
            switch sol.solver
                case 'bvp4c'
                    ocLim=dynprimitive(hatx);
                case 'gbvp4c'
                    ocLim=gdynprimitive(dynprimitive(hatx));
            end
        case 'l'
            lcidx=find(sol.solverinfo.limitcycleindex==sol.solverinfo.octrajectory2limset(ii,3));
            lcleftarcindex=leftarcindex(lcidx);
            lcrightarcindex=rightarcindex(lcidx);
            ocLC.octrajectory.x=sol.x(lcleftarcindex(1):lcrightarcindex(end))-sol.x(lcleftarcindex(1));
            ocLC.octrajectory.y=sol.y(:,lcleftarcindex(1):lcrightarcindex(end));
            ocLC.octrajectory.arcarg=sol.solverinfo.arcarg{hetorder+sol.solverinfo.octrajectory2limset(ii,3)};
            ocLC.octrajectory.arcinterval=sol.solverinfo.arcinterval{hetorder+sol.solverinfo.octrajectory2limset(ii,3)};
            ocLC.octrajectory.timehorizon=ocLC.octrajectory.arcinterval(end);
            ocLC.octrajectory.arcposition=[lcleftarcindex;lcrightarcindex];
            ocLC.octrajectory.arcposition=ocLC.octrajectory.arcposition-ocLC.octrajectory.arcposition(1,1)+1;
            ocLC.octrajectory.modelparameter=sol.modelparameter;
            ocLC.octrajectory.modelname=sol.modelname;
            ocLC.octrajectory.linearization=sol.solverinfo.monodromy;
            ocLC.period=sol.solverinfo.timehorizon(2);
            switch sol.solver
                case 'bvp4c'
                    ocLim=dynprimitive(ocLC);
                case 'gbvp4c'
                    ocLim=gdynprimitive(dynprimitive(ocLC));
            end
    end
    switch sol.solver
        case 'bvp4c'
            ocAsymM{ii}=ocasymptotic(octrajectory(ocTrj),ocLim);
        case 'gbvp4c'
            ocAsymM{ii}=ocgasymptotic(ocgtrajectory(ocTrj,sol.solverinfo.odenum(trjidx)),ocLim);
    end

end
ocAsymM=ocmultipath(ocAsymM);