function ocMP=sol2multipath(sol,ocEP,timehorizon)

if nargin<=2
    timehorizon=sol.timehorizon;
end
midx=numel(timehorizon);
if nargin==1
    ocEP=cell(1,midx);
end
pos=[];
ff=[];
for ii=1:midx
    ff=[ff find(sol.arcinterval==timehorizon(ii))];
%     if isempty(pos)
%         pos=ff(1)-1;
%     else
%         ff(ff<=pos(end)+1)=[];
%         pos=[pos ff(1)-ii];
%     end
end
ff=unique(ff);
pos=[0 cumsum(ff-[1 ff(1:midx-1)+1])];
arcnum=diff(pos);
totalarcidx=[0 cumsum(arcnum)];
lastarcposition=0;
for ii=1:midx
        arcint=sol.arcposition(:,pos(ii)+1:pos(ii+1));
        arcintbd=[arcint(1) arcint(end)];
        msol(ii).x=sol.x(arcintbd(1):arcintbd(2))-sol.x(arcintbd(1));
        msol(ii).y=sol.y(:,arcintbd(1):arcintbd(2));
        if ~isempty(sol.solverinfo)
            if isfield(sol.solverinfo,'inftimetransformation')
                msol(ii).solverinfo.inftimetransformation=sol.solverinfo.inftimetransformation(ii);
            end
            if isfield(sol.solverinfo,'pathtype')
                msol(ii).solverinfo.pathtype=sol.solverinfo.pathtype{ii};
            end
            if isfield(sol.solverinfo,'yp')
                msol(ii).solverinfo.yp=sol.solverinfo.yp(:,arcintbd(1):arcintbd(2));
            end
            msol(ii).solverinfo.multiarccalc=sol.solverinfo.multiarccalc;
        end
        msol(ii).solver=sol.solver;
        msol(ii).arcposition=arcint-lastarcposition;
        msol(ii).arcarg=sol.arcarg(totalarcidx(ii)+1:totalarcidx(ii+1));
        msol(ii).arcinterval=sol.arcinterval(ii+[totalarcidx(ii):totalarcidx(ii+1)]);
        msol(ii).timehorizon=sol.timehorizon(ii);
        if iscell(sol.modelparameter)
            msol(ii).modelparameter=sol.modelparameter{ii};
            msol(ii).modelname=sol.modelname{ii};
        else
            msol(ii).modelparameter=sol.modelparameter;
            msol(ii).modelname=sol.modelname;
        end
        lastarcposition=arcintbd(2);
        if isempty(ocEP{ii})
            ocTrj{ii}=octrajectory(msol(ii));
        else
            ocTrj{ii}=ocasymptotic(octrajectory(msol(ii)),ocEP{ii});
        end
end
ocMP=ocmultipath(ocTrj);