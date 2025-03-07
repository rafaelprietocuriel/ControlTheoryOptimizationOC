function mmTrj2=sol2mmultipath(sol,mmObj,varargin)

if isempty(sol)
    mmTrj2=mmultipath([]);
    return
end
if ~isfield(sol.solverinfo,'numberofstages')
    numberofstages=[];
else
    numberofstages=sol.solverinfo.numberofstages;
end
diffidx=find(diff(sol.x)==0);
leftidx=[1 diffidx+1];
rightidx=[diffidx length(sol.x)];
arcposition=[leftidx;rightidx];
if isfield(sol.solverinfo,'indifferenceorder')
    indifferenceorder=sol.solverinfo.indifferenceorder;
    partstructure=sol.solverinfo.partstructure;
    userinfo.indifferenceorder=indifferenceorder;
    userinfo.partstructure=partstructure;
else
    indifferenceorder=1;
    userinfo=[];
    partstructure=1:length(leftidx);
end

leftidx=[1 find(diff(sol.x)==0)+1];
rightidx=[leftidx(2:end)-1 length(sol.x)];
sol.arcposition=[leftidx;rightidx];
submodname=submodelname(mmObj);
mmTrj=cell(1,numberofstages);
arcpositioncorrect=0;
if indifferenceorder>1
    mmTrj2=cell(1,indifferenceorder);
    for jj=1:indifferenceorder
        ctr=0;
        for ii=partstructure{jj}
            ctr=ctr+1;
            absarcindex=arcposition(1,ii):arcposition(2,ii);
            solStruct.x=sol.x(absarcindex)-sol.x(absarcindex(1));
            solStruct.y=sol.y(:,absarcindex);
            solStruct.arcarg=sol.arcarg(ii);
            solStruct.arcposition=arcposition(:,ii)-arcpositioncorrect;
            arcpositioncorrect=arcposition(2,ii);
            solStruct.arcinterval=sol.arcinterval([ii ii+1]);
            solStruct.modelparameter=sol.modelparameter{ii};
            solStruct.modelname=submodname{ctr};
            solStruct.x0=solStruct.arcinterval(1);
            solStruct.solver=sol.solver;
            switch sol.solver
                case 'bvp4c'
                    solStruct.solverinfo.yp=sol.solverinfo.yp(:,absarcindex);
            end
            solStruct.timehorizon=solStruct.arcinterval(end);
            solStruct.solverinfo.coeff=[];
            solStruct.solverinfo.tangent=[];
            solStruct.solverinfo.tmesh=[];
            solStruct.solverinfo.tmesh=[];
            mmTrj{ctr}=octrajectory(solStruct);
        end
        mmTrj2{jj}=mmultipath(mmTrj);
        sol.arcinterval(1)=[];
    end
else
    for ii=1:numberofstages
        %arc4part=find(sol.solverinfo.solutionindex==ii);
        arc4part=sol.solverinfo.partposition(1,ii):sol.solverinfo.partposition(2,ii);
        absarcindex=(sol.arcposition(1,arc4part(1)):sol.arcposition(2,arc4part(end)));
        solStruct.x=sol.x(absarcindex)-sol.x(absarcindex(1));
        solStruct.y=sol.y(:,absarcindex);
        solStruct.arcarg=sol.arcarg(arc4part);
        solStruct.arcposition=sol.arcposition(:,arc4part)-arcpositioncorrect;
        arcpositioncorrect=sol.arcposition(2,arc4part(end));
        solStruct.arcinterval=sol.arcinterval([arc4part arc4part(end)+1]);
        solStruct.modelparameter=sol.modelparameter{ii};
        solStruct.modelname=submodname{ii};
        solStruct.x0=solStruct.arcinterval(1);
        solStruct.solver=sol.solver;
        switch sol.solver
            case 'bvp4c'
                solStruct.solverinfo.yp=sol.solverinfo.yp(:,absarcindex);
        end
        solStruct.timehorizon=solStruct.arcinterval(end);
        solStruct.solverinfo=sol.solverinfo;
        solStruct.solverinfo.coeff=[];
        solStruct.solverinfo.tangent=[];
        solStruct.solverinfo.tmesh=[];
        solStruct.solverinfo.tmesh=[];
        mmTrj{ii}=octrajectory(solStruct);
    end
    mmTrj2=mmultipath(mmTrj);
end
% if ~isempty(userinfo)
%     mmTrj=mmultipath(mmTrj2,'userinfo',userinfo);
% else
%     mmTrj=mmultipath(mmTrj2);
% end