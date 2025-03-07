function ocMP=sol2ocmultipath(sol,OCMATINDIF,ocEPValue,varargin)

if isocmultipath(sol)
    ocMP=sol;
    return
end

if nargin<=2
    ocEPValue=[];
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
        if isfield(sol.solverinfo,'switchtimecoordinate')
            numarc(ii)=[1+length(sol.solverinfo.switchtimecoordinate{ii})];
        else
            numarc(ii)=[1+length(sol.solverinfo.switchtimecoord{ii})];
        end
        OCMATINDIF.arccoord{ii}=[1:numarc(ii)]+arcoffset;
        arcoffset=arcoffset+numarc(ii);
    end
    OCMATINDIF.cumsumnumarc=cumsum(numarc);
    OCMATINDIF.inftimetransformation=[];%sol.solverinfo.inftimetransformation;
    if isfield(sol.solverinfo,'pathtype')
        OCMATINDIF.pathtype=sol.solverinfo.pathtype;
    end
end

ocAsym=cell(1,OCMATINDIF.indifferenceorder);
cumsumnumarc=[0 OCMATINDIF.cumsumnumarc];
for ii=1:OCMATINDIF.indifferenceorder
    absarcindex=(sol.arcposition(1,cumsumnumarc(ii)+1):sol.arcposition(2,cumsumnumarc(ii+1)));
    solStruct.x=sol.x(absarcindex)-sol.x(absarcindex(1));
    solStruct.y=sol.y(:,absarcindex);
    arccoord=find(sol.solverinfo.solutionindex==ii);
    solStruct.arcarg=sol.arcarg(arccoord);
    solStruct.arcposition=sol.arcposition(:,cumsumnumarc(ii)+1:cumsumnumarc(ii+1))-sol.arcposition(1,cumsumnumarc(ii)+1)+1;
    solStruct.arcinterval=sol.arcinterval(cumsumnumarc(ii)+ii:cumsumnumarc(ii+1)+ii);
    if isfield(sol.solverinfo,'static') && sol.solverinfo.static
        solStruct.solverinfo.static=sol.solverinfo.static;
        solStruct.solverinfo.staticparameterindex=sol.solverinfo.staticparameterindex;
        solStruct.solverinfo.staticparametercoordinate=sol.solverinfo.staticparametercoordinate{ii};
    else
        solStruct.solverinfo.static=0;
    end
    if solStruct.solverinfo.static
        solStruct.modelparameter=sol.modelparameter{ii};
    else
        solStruct.modelparameter=sol.modelparameter;
    end
    solStruct.modelname=sol.modelname;
    solStruct.x0=sol.x0;
    solStruct.solver=sol.solver;
    switch sol.solver
        case 'bvp4c'
            %solStruct.solverinfo.yp=sol.solverinfo.yp(:,absarcindex);
        case 'bvp6c'
            solStruct.solverinfo.yp=sol.solverinfo.yp(:,absarcindex);
            solStruct.solverinfo.ypmid=sol.solverinfo.ypmid(:,absarcindex);
        case 'gbvp4c'
            solStruct.solverinfo.odenum=sol.solverinfo.odenum(arccoord);
            
    end
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
    if isfield(sol.solverinfo,'stateconstraint')
        if sol.solverinfo.stateconstraint{ii}
            solStruct.solverinfo.stateconstraint=sol.solverinfo.stateconstraint{ii};
            solStruct.solverinfo.entryindex=sol.solverinfo.entryindex{ii};
            solStruct.solverinfo.jumpid=sol.solverinfo.jumpid{ii};
            solStruct.solverinfo.jumpargument=sol.solverinfo.parameters(sol.solverinfo.entrytimecoordinate{ii});
        else
            solStruct.solverinfo.stateconstraint=sol.solverinfo.stateconstraint{ii};
            solStruct.solverinfo.entryindex=[];
            solStruct.solverinfo.jumpargument=[];
        end
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
    if isfield(sol.solverinfo,'equilibriumcoordinate')
        switch solStruct.solver
            case 'gbvp4c'
                etype='gdynprimitive';
            otherwise
                etype='dynprimitive';
        end
        ocEPValue{ii}=generateequilibrium(sol.solverinfo.parameters(sol.solverinfo.equilibriumcoordinate{ii}),solStruct.arcarg(end),solStruct.modelparameter,solStruct.modelname,etype);
    end
    if (~iscell(ocEPValue) && isempty(ocEPValue)) || isempty(ocEPValue{ii})
        switch solStruct.solver
            case 'gbvp4c'
                ocAsym{ii}=ocgtrajectory(solStruct);
            otherwise
                ocAsym{ii}=octrajectory(solStruct);
        end
    else
        if isdynprimitive(ocEPValue{ii})
            dynPrim=ocEPValue{ii};
        elseif isa(ocEPValue{ii},'gdynprimitive')
            dynPrim=ocEPValue{ii};
        else
            if isfield(OCMATINDIF,'limitset')
                dynPrim=OCMATINDIF.limitset{ii};
            else
                dynPrim=generateequilibrium(OCMATINDIF.saddlepoint{ii},solStruct.arcarg(1),solStruct.modelparameter,solStruct.modelname);
            end
        end
        switch solStruct.solver
            case 'gbvp4c'
                ocAsym{ii}=ocgasymptotic(ocgtrajectory(solStruct,solStruct.solverinfo.odenum),dynPrim);
            otherwise
                ocAsym{ii}=ocasymptotic(octrajectory(solStruct),dynPrim);
        end
    end
end
ocMP=ocmultipath(ocAsym);