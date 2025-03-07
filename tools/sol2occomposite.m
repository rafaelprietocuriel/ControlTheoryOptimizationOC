function comPath=sol2occomposite(sol,varargin)

diffidx=find(diff(sol.x)==0);
leftidx=[1 diffidx+1];
rightidx=[diffidx length(sol.x)];
arcposition=[leftidx;rightidx];
try
    indifferenceorder=sol.solverinfo.indifferenceorder;
catch
    indifferenceorder=length(sol.modelparameter);
end
arcpositioncorrect=0;
actintervalcounter=0;
comPath=cell(1,indifferenceorder);
for ii=1:indifferenceorder
    try
        indifferenceindex=find(sol.solverinfo.indifferenceindex==ii);
    catch
        indifferenceindex=ii;
    end
    relarcpos=arcposition(:,indifferenceindex);
    try
        partstructure=sol.solverinfo.partindex(indifferenceindex);
    catch
        partstructure=1;
    end
        partnumber=partstructure(end);
    actintervalcounter0=actintervalcounter+1;
    actintervalcounter=actintervalcounter+partnumber+1;
    actarcinterval=sol.arcinterval(actintervalcounter0:actintervalcounter);
    modelparameter=sol.modelparameter(indifferenceindex);
    actarg=sol.arcarg(indifferenceindex);
    ctr=0;
    mmTrj=cell(1,partnumber);
    for jj=1:partnumber
        ctr=ctr+1;
        partindex=find(partstructure==jj);
        partarcpos=relarcpos(:,partindex);
        arcindex=partarcpos(1,1):partarcpos(2,end);
        solStruct.x=sol.x(arcindex)-sol.x(arcindex(1));
        solStruct.y=sol.y(:,arcindex);
        solStruct.arcarg=actarg(partindex);
        solStruct.arcposition=partarcpos-arcpositioncorrect;
        arcpositioncorrect=partarcpos(2,end);
        solStruct.arcinterval=actarcinterval(partindex(1):partindex(end)+1);
        solStruct.timehorizon=solStruct.arcinterval(end);
        solStruct.modelparameter=modelparameter{jj};
        solStruct.solverinfo=sol.solverinfo;
        solStruct.solverinfo=rmfield(solStruct.solverinfo,{'coeff','tangent','tmesh'});
        mmTrj{ctr}=octrajectory(solStruct);
    end
    comPath{ii}=mmultipath(mmTrj);
end
comPath=occomposite(comPath);