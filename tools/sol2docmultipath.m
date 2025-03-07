function docMP=sol2docmultipath(sol,ocFP)

if nargin==1
    ocFP=[];
end
indifferenceorder=size(sol.solverinfo.pathcoord,2);
docAsym=cell(1,indifferenceorder);
numarcidx=[cumsum([1 sol.solverinfo.arcnum(1:end-1)]);cumsum(sol.solverinfo.arcnum)];
totalx=[sol.x0 sol.x];
totaly=[sol.y0 sol.y];
for ii=1:indifferenceorder
    idx=sol.solverinfo.pathcoord(1,ii):sol.solverinfo.pathcoord(2,ii);
    solStruct.x=totalx(idx);
    solStruct.y=totaly(:,idx);
    solStruct.arcarg=sol.arcarg(sol.solverinfo.arcargcoord(1,ii):sol.solverinfo.arcargcoord(2,ii));
    if ii>1
        solStruct.arcposition=sol.arcposition(:,numarcidx(1,ii):numarcidx(2,ii))-(sol.arcposition(2,numarcidx(2,ii-1))+1);
    else
        solStruct.arcposition=sol.arcposition(:,numarcidx(1,ii):numarcidx(2,ii));
    end
    solStruct.modelparameter=sol.modelparameter;
    solStruct.modelname=sol.modelname;
    solStruct.x0=solStruct.x(1);
    solStruct.x(1)=[];
    solStruct.y0=solStruct.y(:,1);
    solStruct.y(:,1)=[];

    solStruct.timehorizon=sol.timehorizon(ii);
    solStruct.solverinfo.pathtype=sol.solverinfo.pathtype{ii};
    solStruct.solverinfo.coeff=[];
    solStruct.solverinfo.tangent=[];
    solStruct.solverinfo.tmesh=[];
    if ~isempty(ocFP)
        docAsym{ii}=docasymptotic(doctrajectory(solStruct),ocFP{ii});
    else
        docAsym{ii}=doctrajectory(solStruct);
    end
end
docMP=docmultipath(docAsym);