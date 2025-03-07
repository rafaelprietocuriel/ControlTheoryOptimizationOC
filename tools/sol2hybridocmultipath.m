function ocMP=sol2hybridocmultipath(sol,IOCMATINDIF)

if nargin==1 || isempty(IOCMATINDIF)
    IOCMATINDIF.indifferenceorder=length(sol.solverinfo.pathtype);
    numarc=zeros(1,IOCMATINDIF.indifferenceorder);
    arcoffset=0;
    for ii=1:IOCMATINDIF.indifferenceorder
        numarc(ii)=[1+length(sol.solverinfo.switchtimecoord{ii})];
        IOCMATINDIF.arccoord{ii}=[1:numarc(ii)]+arcoffset;
        arcoffset=numarc(ii);
    end
    IOCMATINDIF.cumsumnumarc=cumsum(numarc);
    IOCMATINDIF.inftimetransformation=[];%sol.solverinfo.inftimetransformation;
    IOCMATINDIF.pathtype=sol.solverinfo.pathtype;
    IOCMATINDIF.initialdepvarcoord=sol.solverinfo.initialdepvarcoord;
    IOCMATINDIF.enddepvarcoord=sol.solverinfo.enddepvarcoord;
    IOCMATINDIF.jumpcoord=sol.solverinfo.jumpcoord;
end
hocTrj=cell(1,IOCMATINDIF.indifferenceorder);
cumsumnumarc=[0 IOCMATINDIF.cumsumnumarc];
for ii=1:IOCMATINDIF.indifferenceorder
    X0=sol.parameters(IOCMATINDIF.initialdepvarcoord{ii});
    XT=sol.parameters(IOCMATINDIF.enddepvarcoord{ii});
    absarcindex=[sol.arcposition(1,cumsumnumarc(ii)+1):sol.arcposition(2,cumsumnumarc(ii+1))];
    solStruct.x=sol.x(absarcindex)-sol.x(absarcindex(1));
    solStruct.x=[solStruct.x(1) solStruct.x solStruct.x(end)];
    if length(X0)<size(sol.y,1)
        solStruct.y=[[X0;0] sol.y(:,absarcindex) [XT;sol.y(end,absarcindex(end))]];
    else
        solStruct.y=[X0 sol.y(:,absarcindex) XT];
    end
    solStruct.arcarg=sol.arcarg(IOCMATINDIF.arccoord{ii});
    solStruct.jumparg=sol.solverinfo.jumparg(IOCMATINDIF.jumpcoord{ii});
    arcposition=find(diff(solStruct.x(2:end-1))==0)+1;
    solStruct.arcposition=[1 arcposition+1;arcposition length(solStruct.x)];
    solStruct.arcposition(1,1)=solStruct.arcposition(1,1)+1;
    solStruct.arcposition(2,end)=solStruct.arcposition(2,end)-1;
    solStruct.arcinterval=sol.arcinterval(cumsumnumarc(ii)+ii:cumsumnumarc(ii+1)+ii);
    solStruct.modelparameter=sol.modelparameter;
    solStruct.modelname=sol.modelname;
    solStruct.x0=sol.x0;
    solStruct.solver=sol.solver;
    switch sol.solver
        case 'bvp4c'
            solStruct.solverinfo.yp=sol.solverinfo.yp(:,absarcindex);
        case 'bvp6c'
            solStruct.solverinfo.yp=sol.solverinfo.yp(:,absarcindex);
            solStruct.solverinfo.ypmid=sol.solverinfo.ypmid(:,absarcindex);
    end
    if isfield(sol.solverinfo,'switchtimecoord')
        solStruct.solverinfo.switchtimecoord=sol.solverinfo.switchtimecoord{ii};
    end
    solStruct.timehorizon=sol.timehorizon(ii);
    solStruct.solverinfo.pathtype=IOCMATINDIF.pathtype{ii};
    solStruct.solverinfo.inftimetransformation=[];
    solStruct.solverinfo.coeff=[];
    solStruct.solverinfo.tangent=[];
    solStruct.solverinfo.tmesh=[];
    hocTrj{ii}=hybridoctrajectory(solStruct);
end
ocMP=hybridocmultipath(hocTrj);