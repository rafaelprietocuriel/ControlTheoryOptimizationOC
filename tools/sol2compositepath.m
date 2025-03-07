function ocMP=sol2compositepath(sol,OCMATINDIF,varargin)

if isocmultipath(sol)
    ocMP=sol;
    return
end

if nargin==1 || isempty(OCMATINDIF)
    if isfield(sol.solverinfo,'indifferenceorder')
        OCMATINDIF.indifferenceorder=sol.solverinfo.indifferenceorder;
    else
        OCMATINDIF.indifferenceorder=length(sol.solverinfo.pathtype);
    end
    numarc=zeros(1,OCMATINDIF.indifferenceorder);
    arcoffset=0;
    for ii=1:OCMATINDIF.indifferenceorder
        numarc(ii)=[1+length(sol.solverinfo.switchtimecoord{ii})];
        OCMATINDIF.arccoord{ii}=[1:numarc(ii)]+arcoffset;
        arcoffset=arcoffset+numarc(ii);
    end
    OCMATINDIF.cumsumnumarc=cumsum(numarc);
    OCMATINDIF.inftimetransformation=[];%sol.solverinfo.inftimetransformation;
    if isfield(sol.solverinfo,'pathtype')
        OCMATINDIF.pathtype=sol.solverinfo.pathtype;
    end
end
ocTrj=cell(1,OCMATINDIF.indifferenceorder);
cumsumnumarc=[0 OCMATINDIF.cumsumnumarc];
for ii=1:OCMATINDIF.indifferenceorder
    absarcindex=(sol.arcposition(1,cumsumnumarc(ii)+1):sol.arcposition(2,cumsumnumarc(ii+1)));
    solStruct.x=sol.x(absarcindex)-sol.x(absarcindex(1));
    solStruct.y=sol.y(:,absarcindex);
    solStruct.arcarg=sol.arcarg(OCMATINDIF.arccoord{ii});
    solStruct.arcposition=sol.arcposition(:,cumsumnumarc(ii)+1:cumsumnumarc(ii+1))-sol.arcposition(1,cumsumnumarc(ii)+1)+1;
    solStruct.arcinterval=sol.arcinterval(cumsumnumarc(ii)+ii:cumsumnumarc(ii+1)+ii);
    solStruct.modelparameter=sol.modelparameter{ii};
    solStruct.modelname=sol.modelname{ii};
    solStruct.x0=sol.x0;
    solStruct.solver=sol.solver;
    if isfield(sol.solverinfo,'switchtimecoord')
        solStruct.solverinfo.switchtimecoord=sol.solverinfo.switchtimecoord{ii};
    end
    solStruct.timehorizon=sol.timehorizon(ii);
    try
        solStruct.solverinfo.continuationparameter=sol.solverinfo.continuationparameter;
    catch
        solStruct.solverinfo.continuationparameter=[];
    end
    if isfield(OCMATINDIF,'pathtype')
        solStruct.solverinfo.pathtype=OCMATINDIF.pathtype{ii};
    end
    if isfield(sol.solverinfo,'vjumpargument')
        solStruct.solverinfo.vjumpargument=sol.solverinfo.vjumpargument{ii};
    else
        solStruct.solverinfo.vjumpargument=[];
    end
    if isfield(sol.solverinfo,'vfreetime')
        solStruct.solverinfo.vfreetime=sol.solverinfo.vfreetime{ii};
    else
        solStruct.solverinfo.vfreetime=[];
    end
    solStruct.solverinfo.indifferenceorder=OCMATINDIF.indifferenceorder;
    solStruct.solverinfo.inftimetransformation=[];
    solStruct.solverinfo.coeff=[];
    solStruct.solverinfo.tangent=[];
    solStruct.solverinfo.tmesh=[];
        ocTrj{ii}=octrajectory(solStruct);
end
ocMP=occomposite(ocTrj);