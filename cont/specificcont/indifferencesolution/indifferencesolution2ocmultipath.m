function ocAsymM=indifferencesolution2ocmultipath(ocObj,sol)


switch sol.solverinfo.conttype
    case 'indifferencesolutionp'
        par=sol.modelparameter;
        ocObj=changeparametervalue(ocObj,par);
        arcposition=find(diff(sol.x)==0);
        leftarcindex=[1 arcposition+1];
        rightarcindex=[arcposition numel(sol.x)];

        %
        % for ii=1:sol.solverinfo.equilibriumcounter
        %     hatx.y=sol.solverinfo.parameters(sol.solverinfo.equilibriumcoordinate{ii});
        %     hatx.x=0;
        %     hatx.arcarg=sol.solverinfo.arcarg{3}(end);
        %     hatx.linearization=jacobian(ocObj,hatx);
        %     ocEP=dynprimitive(hatx);
        % end

        indifforder=length(sol.solverinfo.pathtype);
        ocAsymM=cell(1,indifforder);
        arcidx=0;
        for ii=1:indifforder
            trjidx=find(sol.solverinfo.solutionindex==ii);
            trjleftarcindex=leftarcindex(trjidx);
            trjrightarcindex=rightarcindex(trjidx);

            numarc=length(sol.solverinfo.arcarg{ii});
            ocTrj.x=sol.x(trjleftarcindex(1):trjrightarcindex(end))-sol.x(trjleftarcindex(1));
            ocTrj.y=sol.y(:,trjleftarcindex(1):trjrightarcindex(end));
            ocTrj.arcarg=sol.solverinfo.arcarg{ii};
            ocTrj.arcinterval=sol.arcinterval(arcidx+(1:numarc));
            arcidx=arcidx+numarc;
            ocTrj.timehorizon=sol.timehorizon(ii);
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
                    ocAsymM{ii}=ocasymptotic(octrajectory(ocTrj),dynprimitive(hatx));
                case 'l'
                    lcidx=find(sol.solverinfo.limitcycleindex==sol.solverinfo.octrajectory2limset(ii,3));
                    lcleftarcindex=leftarcindex(lcidx);
                    lcrightarcindex=rightarcindex(lcidx);
                    ocLC.octrajectory.x=sol.x(lcleftarcindex(1):lcrightarcindex(end))-sol.x(lcleftarcindex(1));
                    ocLC.octrajectory.y=sol.y(:,lcleftarcindex(1):lcrightarcindex(end));
                    ocLC.octrajectory.arcarg=sol.solverinfo.arcarg{indifforder+sol.solverinfo.octrajectory2limset(ii,3)};
                    ocLC.octrajectory.arcinterval=sol.solverinfo.arcinterval{indifforder+sol.solverinfo.octrajectory2limset(ii,3)};
                    ocLC.octrajectory.timehorizon=ocLC.octrajectory.arcinterval(end);
                    ocLC.octrajectory.arcposition=[lcleftarcindex;lcrightarcindex];
                    ocLC.octrajectory.arcposition=ocLC.octrajectory.arcposition-ocLC.octrajectory.arcposition(1,1)+1;
                    ocLC.octrajectory.modelparameter=sol.modelparameter;
                    ocLC.octrajectory.modelname=sol.modelname;
                    ocLC.octrajectory.linearization=sol.solverinfo.monodromy;
                    ocLC.period=sol.solverinfo.timehorizon(2);

                    ocAsymM{ii}=ocasymptotic(octrajectory(ocTrj),dynprimitive(ocLC));
            end

        end
        ocAsymM=ocmultipath(ocAsymM);
    case 'indifferencesolution2lc'
        par=sol.modelparameter;
        ocObj=changeparametervalue(ocObj,par);
        arcposition=find(diff(sol.x)==0);
        leftarcindex=[1 arcposition+1];
        rightarcindex=[arcposition numel(sol.x)];

        %
        % for ii=1:sol.solverinfo.equilibriumcounter
        %     hatx.y=sol.solverinfo.parameters(sol.solverinfo.equilibriumcoordinate{ii});
        %     hatx.x=0;
        %     hatx.arcarg=sol.solverinfo.arcarg{3}(end);
        %     hatx.linearization=jacobian(ocObj,hatx);
        %     ocEP=dynprimitive(hatx);
        % end

        indifforder=length(sol.solverinfo.pathtype);
        ocAsymM=cell(1,indifforder);
        arcidx=0;
        for ii=1:indifforder
            if ~isempty(sol.solverinfo.pathtype{ii})
                trjidx=find(sol.solverinfo.solutionindex==ii);
                trjleftarcindex=leftarcindex(trjidx);
                trjrightarcindex=rightarcindex(trjidx);

                numarc=length(sol.solverinfo.arcarg{ii});
                ocTrj.x=sol.x(trjleftarcindex(1):trjrightarcindex(end))-sol.x(trjleftarcindex(1));
                ocTrj.y=sol.y(:,trjleftarcindex(1):trjrightarcindex(end));
                ocTrj.arcarg=sol.solverinfo.arcarg{ii};
                ocTrj.arcinterval=sol.arcinterval(arcidx+(1:numarc+1));
                arcidx=arcidx+numarc+1;
                ocTrj.timehorizon=sol.timehorizon(ii);
                ocTrj.arcposition=[trjleftarcindex;trjrightarcindex];
                ocTrj.arcposition=ocTrj.arcposition-ocTrj.arcposition(1,1)+1;
                ocTrj.modelparameter=sol.modelparameter;
                ocTrj.modelname=sol.modelname;
                ocTrj.solverinfo.pathtype=sol.solverinfo.pathtype{ii};
                switch sol.solverinfo.limitsettype{sol.solverinfo.octrajectory2limset(ii,2)}
                    case 'e'
                        hatx.y=sol.solverinfo.parameters(sol.solverinfo.equilibriumcoordinate{sol.solverinfo.octrajectory2limset(ii,2)});
                        hatx.x=0;
                        hatx.arcarg=sol.solverinfo.arcarg{ii}(end);
                        hatx.linearization=jacobian(ocObj,hatx);
                        ocAsymM{ii}=ocasymptotic(octrajectory(ocTrj),dynprimitive(hatx));
                    case 'l'
                        lcidx=find(sol.solverinfo.limitcycleindex==sol.solverinfo.octrajectory2limset(ii,2));
                        lcleftarcindex=leftarcindex(lcidx);
                        lcrightarcindex=rightarcindex(lcidx);
                        ocLC.octrajectory.x=sol.x(lcleftarcindex(1):lcrightarcindex(end))-sol.x(lcleftarcindex(1));
                        ocLC.octrajectory.y=sol.y(:,lcleftarcindex(1):lcrightarcindex(end));
                        ocLC.octrajectory.arcarg=sol.solverinfo.arcarg{indifforder+sol.solverinfo.octrajectory2limset(ii,2)};
                        ocLC.octrajectory.arcinterval=sol.solverinfo.arcinterval{indifforder+sol.solverinfo.octrajectory2limset(ii,2)};
                        ocLC.octrajectory.timehorizon=ocLC.octrajectory.arcinterval(end);
                        ocLC.octrajectory.arcposition=[lcleftarcindex;lcrightarcindex];
                        ocLC.octrajectory.arcposition=ocLC.octrajectory.arcposition-ocLC.octrajectory.arcposition(1,1)+1;
                        ocLC.octrajectory.modelparameter=sol.modelparameter;
                        ocLC.octrajectory.modelname=sol.modelname;
                        ocLC.octrajectory.linearization=sol.solverinfo.monodromy;
                        ocLC.period=sol.solverinfo.timehorizon(2);

                        ocAsymM{ii}=ocasymptotic(octrajectory(ocTrj),dynprimitive(ocLC));
                end
            else
                trjidx=find(sol.solverinfo.solutionindex==ii);
                trjleftarcindex=leftarcindex(trjidx);
                trjrightarcindex=rightarcindex(trjidx);

                numarc=length(sol.solverinfo.arcarg{ii});
                ocLC.octrajectory.x=sol.x(trjleftarcindex(1):trjrightarcindex(end))-sol.x(trjleftarcindex(1));
                ocLC.octrajectory.y=sol.y(:,trjleftarcindex(1):trjrightarcindex(end));
                ocLC.octrajectory.arcarg=sol.solverinfo.arcarg{ii};
                ocLC.octrajectory.arcinterval=sol.solverinfo.arcinterval{ii};
                arcidx=arcidx+numarc;
                ocLC.octrajectory.timehorizon=ocLC.octrajectory.arcinterval(end);
                ocLC.octrajectory.arcposition=[trjleftarcindex;trjrightarcindex];
                ocLC.octrajectory.arcposition=ocLC.octrajectory.arcposition-ocLC.octrajectory.arcposition(1,1)+1;
                ocLC.octrajectory.modelparameter=sol.modelparameter;
                ocLC.octrajectory.modelname=sol.modelname;
                ocLC.period=ocLC.octrajectory.timehorizon;
                ocAsymM{ii}=dynprimitive(ocLC);
            end
        end
        ocAsymM=ocmultipath(ocAsymM);
end