function sol=initocmat_MM_FTE_IS(mmObj,ocTrjMP,continuationtype,continuationindex,continuationtarget,varargin)
%
% initocmat_MM_FTE_IS initialization for indifference solution continuation
% of multimodel solution
%
% SOL=initocmat_MM_FTE_IS(MMOBJ,OCTRJM,CONTINUATIONTYPE,CONTINUATIONINDEX,CONTINUATIONTARGET) the
% continuation process is initialized (started) at a solution OCTRJM (cell
% of mmultipath solutions
%
% We assume that the different solutions are of the same structure of parts
% (stages) and the according models are the same for each corresponding
% part
% EXAMPLE: three parts, indifference order two: if the first solution is a
% solution of the models m_1/par_1, m_2/par_2, m_3/par_3, then the second
% solution is of the same sturcture, same models with same parameter values
%
% We will also allow one part/stage solutions of any indifference order:
% this means that for every solution we allow a different model and or
% different parameter values
%
% FREEPARAMETERVECTOR: since the parameter values or/an its location can be
% different for every part in every indifference solution, we have to
% specify the indices of free parameters for each parametervector.
% Example: two parts, indifference order two: alpha,beta appear in all
% parts at the same location, say 1,2 then the freeparametervector is given
% by a cellarray of length four: {[1 2]},{[1 2]},{[1 2]},{[1 2]} the target
% value is of the same structure
%
%two parts, indifference order two: alpha,beta appear in different parts at
%a different location´, say part1: [1 2], part2: [2 1], then the
%freeparametervector is given by a cellaray of length four: {[1 2]},{[2
%1]},{[1 2]},{[2 1]}
%
%two parts, indifference order two: in the first part we want to continue
%along alpha, beta, [1 2] in the second part along sigma and gamma [3 4]
%then the freeparametervector is given by a cellarray of length four:  
% {[1 2]},{[3 4]},{[1 2]},{[3 4]} the target 
% value is of the same structure
% For a more detailed discussion see INITOCMAT_AE_EP
%
clear global OCMATCONT OCMATINDIF
global OCMATCONT OCMATINDIF
sol=[];

fixendstate=[];
fixinitialstate=[];
freeparametervector=[];
freestatevector=[];
freetimevector=[];

objectivevaluecalc=[];

if isempty(mmObj)
    ocmatmsg('oc model is empty.')
    return
end

% contiuationtype: time, state, parameter
% continuationindex: (time) number of switching/endtime, (state)
% coordinate, (parameter) index or name of parameter
% continuationtarget: value
fixinitialstateidx=find(strcmpi(varargin,'fixinitialstate'));
fixendstateidx=find(strcmpi(varargin,'fixendstate'));
freeparametervectoridx=find(strcmpi(varargin,'freeparametervector'));
freestatevectoridx=find(strcmpi(varargin,'freestatevector'));
freetimevectoridx=find(strcmpi(varargin,'freetimevector'));
objectivevaluecalcidx=find(strcmpi(varargin,'objectivevaluecalc'));
if ~isempty(fixinitialstateidx)
    fixinitialstate=varargin{fixinitialstateidx+1};
end
if ~isempty(fixendstateidx)
    fixendstate=varargin{fixendstateidx+1};
end
if ~isempty(freeparametervectoridx)
    freeparametervector=varargin{freeparametervectoridx+1};
end
if ~isempty(freestatevectoridx)
    freestatevector=varargin{freestatevectoridx+1};
end
if ~isempty(freetimevectoridx)
    freetimevector=varargin{freetimevectoridx+1};
end
if ~isempty(objectivevaluecalcidx)
    objectivevaluecalc=varargin{objectivevaluecalcidx+1};
end

if isempty(objectivevaluecalc)
    objectivevaluecalc=0;
end

if iscell(ocTrjMP)
    userinfo.multiplicity=length(ocTrjMP);
    ctr=0;
    for ii=1:length(ocTrjMP)
        if ismmultipath(ocTrjMP{ii})
            ctr0=ctr+1;
            ctr=ctr+numberofparts(ocTrjMP{ii});
            userinfo.partstructure{ii}=ctr0:ctr;
        else
        end
    end
    ocTrjMP=mmultipath(ocTrjMP,'userinfo',userinfo);
else
end
nummod=numberofmodels(mmObj);
numparts=numberofparts(ocTrjMP);
indifforder=multiplicity(ocTrjMP);
if isempty(indifforder)
    return
end
partstruct=partstructure(ocTrjMP);
if isempty(partstruct)
    return
end
arcstruct=cell(1,numparts);
arcarg=cell(1,numparts);
arcpos=cell(1,numparts);
ctr=0;
for ii=1:numparts
    ctr0=ctr+1;
    arcarg{ii}=arcargument(ocTrjMP(ii));
    arcpos{ii}=arcposition(ocTrjMP(ii));
    ctr=ctr+length(arcarg{ii});
    arcstruct{ii}=ctr0:ctr;
end
OCMATINDIF.partstructure=partstruct;
OCMATINDIF.arcstructure=arcstruct;
OCMATINDIF.arcargument=arcarg;
OCMATINDIF.arcposition=arcpos;
OCMATINDIF.numberofparts=numparts;
OCMATINDIF.indifferenceorder=indifforder;
OCMATINDIF.optimalswitchingtime=optimaltransition(ocTrjMP);
OCMATINDIF.objectivevaluecalc=objectivevaluecalc;
% initialize global variable (OCMATINDIF) for specific continuation process 
% initialize functions 
% empty function handles mean either that they have to be computed
% numerically (derivatives, jacobian, hessian) or that standard functions are
% used (e.g. asymptotic boundary condition, initial condition)
% mode and path specific variables
if numparts==nummod
    OCMATINDIF.parametervalue=parametervalue(mmObj);
elseif numparts==nummod*indifforder
    OCMATINDIF.parametervalue=repmat(parametervalue(mmObj),1,indifforder);
end

% freestatevector is a matrix, where each column is a free vector and the
% number of rows is equal to the number of states
OCMATINDIF.freestatevector=freestatevector;
% freeparametervector: we have the same assumptions as for the parameter
% continuation described above, thus for each parametervector a cellarray
% of such a structure has to be provided.
% EXAMPLE: {{'alpha,beta','alpha,beta','alpha,beta'},{'alpha,sigma','alpha,sigma','alpha,sigma'}}
%

OCMATINDIF.freeparametervector=freeparametervector;
if ~isempty(freeparametervector)
    for ii=1:length(freeparametervector)
        if iscell(freeparametervector{ii}) && length(freeparametervector{ii})==nummod
            freeparametervector{ii}=repmat(freeparametervector{ii},1,indifforder);
        elseif ischar(freeparametervector{ii})
            freeparametervector{ii}=repmat({freeparametervector{ii}},1,numparts);
        end
        if numparts==nummod
            OCMATINDIF.freeparametervectorcoordinate=parameterindex(mmObj,freeparametervector{ii}(partstruct));
        elseif numparts==nummod*indifforder
            ctr=0;
            ctr2=0;
            for jj=1:indifforder
                ctr20=ctr2+1;
                ctr2=ctr2+nummod;
                % these are the free indices for each model parameter set
                actparindex=parameterindex(mmObj,freeparametervector{ii}(partstruct{jj}));
                OCMATINDIF.freeparametervectorcoordinate{ii}(:,[ctr20:ctr2])=reshape([actparindex{:}],[],nummod);
                % we assume that e.g. [alpha, beta]
                OCMATINDIF.freeparametervector{ii}(:,[ctr20:ctr2])=ones(size(OCMATINDIF.freeparametervectorcoordinate{ii}(:,[ctr20:ctr2])));
                for kk=1:nummod
                    ctr=ctr+1;
                    OCMATINDIF.freeinitialparametervalue{ii}(:,ctr)=OCMATINDIF.parametervalue{ctr}(OCMATINDIF.freeparametervectorcoordinate{ii}(:,ctr20+kk-1));
                end
            end
        end
    end
end
OCMATINDIF.freetimevector=freetimevector;
if ~isempty(freetimevector)
end

% initialize global variable (OCMATCONT) for general continuation process
OCMATCONT.modelname=modelname(mmObj);
OCMATCONT.modelfunc=modelspecificfunc(mmObj,'4IndifferenceSolutionContinuation');

ctr=0;
for ii=1:nummod
    if numparts==nummod
    elseif numparts==nummod*indifforder
        for jj=1:indifforder
            ctr=ctr+1;
            funch{ctr}=OCMATCONT.modelfunc{ii}(); % model specific function handles for saddle path continuation
        end
    end
end
ctr=0;
for ii=1:numparts
    OCMATINDIF.canonicalsystem{ii}=funch{ii}{1};
    OCMATINDIF.canonicalsystemjacobian{ii}=funch{ii}{2}{1};
    OCMATINDIF.canonicalsystemparameterjacobian{ii}=funch{ii}{2}{2};
    OCMATINDIF.canonicalsystemhessian{ii}=funch{ii}{3}{1};
    OCMATINDIF.canonicalsystemparameterhessian{ii}=funch{ii}{3}{2};

    % function for the boundary conditions
    OCMATINDIF.bcinitial{ii}=funch{ii}{5}{1};
    OCMATINDIF.bcasymptotic{ii}=funch{ii}{5}{2};
    OCMATINDIF.bctransversality{ii}=funch{ii}{5}{3};
    OCMATINDIF.bcindifference{ii}=funch{ii}{5}{5};
    OCMATINDIF.salvagevalue{ii}=funch{ii}{5}{6};
    OCMATINDIF.bcconnectingparts{ii}=funch{ii}{5}{7};
    OCMATINDIF.bcoptimalconnectingparts{ii}=funch{ii}{5}{8};
    
    % function for Jacobian
    OCMATINDIF.bcjacobianinitial{ii}=funch{ii}{6}{1};
    OCMATINDIF.bcjacobiantransversality{ii}=funch{ii}{6}{2};

    % function describing the hybrid structure of the problem
    OCMATINDIF.hybridinfo{ii}=funch{ii}{7}{1};
    OCMATINDIF.domain{ii}=funch{ii}{7}{2};
    OCMATINDIF.guard{ii}=funch{ii}{7}{3};
    OCMATINDIF.reset{ii}=funch{ii}{7}{4};
    OCMATINDIF.switchtime{ii}=funch{ii}{7}{5};
    OCMATINDIF.jacobianguard{ii}=funch{ii}{7}{7};
    OCMATINDIF.jacobianreset{ii}=funch{ii}{7}{8};
    OCMATINDIF.domaindiscretization{ii}=funch{ii}{7}{9};

    OCMATINDIF.objectivefunction{ii}=funch{ii}{8}{1};
    OCMATINDIF.objectivefunctionjacobian{ii}=funch{ii}{8}{2};
    OCMATINDIF.objectivefunctionparameterjacobian{ii}=funch{ii}{8}{3};
    OCMATINDIF.objectivefunctionderivativetime{ii}=funch{ii}{8}{4};
    
    if numparts==nummod
        OCMATINDIF.autonomous(ii)=isautonomous(mmObj(ii));
        if ~isautonomous(mmObj(ii))
            OCMATINDIF.canonicalsystemderivativetime{ii}=funch{ii}{2}{3};
        end
    elseif numparts==nummod*indifforder
        ctr=ctr+1;
        if ctr>nummod
            ctr=ctr-nummod;
        end
        OCMATINDIF.autonomous(ii)=isautonomous(mmObj(ctr));
        if ~isautonomous(mmObj(ctr))
            OCMATINDIF.canonicalsystemderivativetime{ii}=funch{ii}{2}{3};
        end
    end
    % general function
    OCMATINDIF.plotcontinuation{ii}=funch{ii}{11};
    OCMATINDIF.testadmissibility{ii}=funch{ii}{12};
    OCMATINDIF.datapath{ii}=funch{ii}{20};
    OCMATINDIF.saveintermediatefiles{ii}=funch{ii}{21};

    hybridinfo{ii}=OCMATINDIF.hybridinfo{ii}();
    
    for jj=1:numel(hybridinfo{ii}.arcarg)
        domaindata{ii}(jj)=OCMATINDIF.domain{ii}(hybridinfo{ii}.arcarg(jj));
    end
    for jj=1:numel(domaindata{ii})
        OCMATCONT.DOMAINDDATA{ii}(jj).numode=domaindata{ii}(jj).odedim;
        OCMATCONT.DOMAINDDATA{ii}(jj).numae=domaindata{ii}(jj).aedim;
        OCMATCONT.DOMAINDDATA{ii}(jj).daeorder=domaindata{ii}(jj).daeorder;
        OCMATCONT.DOMAINDDATA{ii}(jj).numeq=numel(domaindata{ii}(jj).daeorder);%number of equations
        if objectivevaluecalc
            OCMATCONT.DOMAINDDATA{ii}(jj).numode=domaindata{ii}(1).odedim+1;
            if jj==1
                OCMATCONT.DOMAINDDATA{ii}(jj).numeq=OCMATCONT.DOMAINDDATA{ii}(1).numeq+1;%number of equations
            else
                OCMATCONT.DOMAINDDATA{ii}(jj).numeq=OCMATCONT.DOMAINDDATA{ii}(1).numeq;%number of equations
            end
        end
        OCMATCONT.DOMAINDDATA{ii}(jj).eqcoord=1:OCMATCONT.DOMAINDDATA{ii}(jj).numeq;
        OCMATCONT.DOMAINDDATA{ii}(jj).odecoord=1:domaindata{ii}(jj).odedim;
        OCMATCONT.DOMAINDDATA{ii}(jj).aecoord=domaindata{ii}(jj).odedim+(1:domaindata{ii}(jj).aedim);
    end
end

OCMATINDIF.initialstate=state(mmObj(1),ocTrjMP(1),1);
OCMATINDIF.initialstate=OCMATINDIF.initialstate(:,1);
OCMATINDIF.fixinitialstate=fixinitialstate;
if ~isempty(fixinitialstate)
    OCMATINDIF.fixinitialstatecoordinate=fixinitialstate;
end
OCMATINDIF.fixendstate=fixendstate;
arcint=arcinterval(ocTrjMP);
arcint=[arcint{:}];
arcint(diff(arcint)==0)=[];
OCMATINDIF.initialarcinterval=arcint;
OCMATINDIF.continuationindex=continuationindex;

switch continuationtype
    case 'initialstate'
        OCMATINDIF.continuationvector=continuationtarget(:)-OCMATINDIF.initialstate(continuationindex);
        continuationparameter=0;
        OCMATINDIF.continuationtype=0;
    case 'endstate'
        OCMATINDIF.continuationvector=continuationtarget(:)-OCMATINDIF.endstate(continuationindex);
        continuationparameter=0;
        OCMATINDIF.continuationtype=0;

    case {'switchingtime','time'}
        OCMATINDIF.continuationvector=continuationtarget(:).'-OCMATINDIF.initialarcinterval(continuationindex);
        continuationparameter=0;
        OCMATINDIF.continuationindex=continuationindex;
        OCMATINDIF.continuationtype=1;
    case 'parameter'
        OCMATINDIF.continuationtype=2;
        parindex=parameterindex(mmObj,continuationindex);
        if isnumeric(continuationtarget)
            if length(continuationtarget)==nummod
                continuationtarget=num2cell(repmat(continuationtarget(:).',1,indifforder));
            else
                continuationtarget=repmat({continuationtarget},1,numparts);
            end
        end
        if iscell(continuationtarget) && length(continuationtarget)==1

            OCMATINDIF.continuationvector=cell(1,numparts*indifforder);
            ctr=0;
            for ii=1:indifforder
                for jj=1:numparts
                    ctr=ctr+1;
                    OCMATINDIF.continuationvector{ctr}=continuationtarget{jj}-OCMATINDIF.parametervalue{jj}(parindex{jj});
                    OCMATINDIF.initialparameter{ctr}=OCMATINDIF.parametervalue{jj}(parindex{jj});
                    OCMATINDIF.parameterindex{ctr}=parindex{jj};
                end
            end
        else
            OCMATINDIF.continuationvector=cell(1,nummod*indifforder);
            ctr=0;
            for ii=1:indifforder
                for jj=1:nummod
                    ctr=ctr+1;
                    OCMATINDIF.continuationvector{ctr}=continuationtarget{ctr}-OCMATINDIF.parametervalue{jj}(parindex{jj});
                    OCMATINDIF.initialparameter{ctr}=OCMATINDIF.parametervalue{jj}(parindex{jj});
                    OCMATINDIF.parameterindex{ctr}=parindex{jj};
                end
            end

        end
        continuationparameter=0;
        OCMATINDIF.continuationindex=parindex{1};
end
scoord=statecoord(mmObj);
cscoord=costatecoord(mmObj);
if objectivevaluecalc
    objpath=objectivevaluepath(ocTrjMP);
    if isempty(objpath)
        ctr=0;
        objpath=cell(1,numparts);
        for ii=1:indifforder
            Oval=0;
            for jj=1:nummod
                ctr=ctr+1;
                if numparts==nummod
                    actmmObj=mmObj(ctr);
                elseif numparts==nummod*indifforder
                    actmmObj=mmObj(jj);
                end
                if length(ocTrjMP(ctr).y(:,1))==2*statenum(actmmObj)
                    o=objectivefunction(actmmObj,ocTrjMP(ctr),1);
                    objpath{ctr}=Oval+[0 cumsum((o(1:end-1)+o(2:end))/2.*diff(time(actmmObj,ocTrjMP(ii),1)))];
                    Oval=ocTrjStruct.y(end,end);
                else
                    objpath{ctr}=ocTrjMP(ctr).y(end,:);
                    Oval=objpath{ctr}(end);
                end
            end
        end
    end
    
    OCMATINDIF.objectivevaluecoord=cscoord{1}(end)+1;
end
ocTrjMP=addobjectivevaluepath(ocTrjMP,objpath,cscoord{1}(end)+1);
ctr=0;
for ii=1:numparts
    if numparts==nummod
        OCMATINDIF.statecostatecoord{ii}=[scoord{ii}(:).' cscoord{ii}(:).'];
    elseif numparts==nummod*indifforder
        ctr=ctr+1;
        if ctr>nummod
            ctr=ctr-nummod;
        end
        OCMATINDIF.statecostatecoord{ii}=[scoord{ctr}(:).' cscoord{ctr}(:).'];
    end
end


sol=generatesolstruct(ocTrjMP);

depvarFirst=zeros(cscoord{1}(end),indifforder);
depvarLast=zeros(cscoord{1}(end),indifforder);
for ii=1:indifforder
    depvar=dependentvar(ocTrjMP(partstruct{ii}(1)));
    depvarFirst(:,ii)=depvar(OCMATINDIF.statecostatecoord{partstruct{ii}(1)},1);
    depvar=dependentvar(ocTrjMP(partstruct{ii}(end)));
    depvarLast(:,ii)=depvar(OCMATINDIF.statecostatecoord{partstruct{ii}(end)},end);
end
OCMATINDIF.initialstate=depvarFirst(scoord{1},:);
OCMATINDIF.endstate=depvarLast(scoord{nummod},:);

optflag=0;
ctrp=0;
for ii=1:indifforder
    switchtimeindex{ii}=[]; % switches between arcs of different constraint number
    parttimeindex{ii}=[]; % switches between different parts/stages
    switchparttimeindex{ii}=[]; % switches between different parts/stages
    optimalparttimeindex{ii}=[]; % switches between different parts/stages
    partcounter=0;
    ctr=0;
    for jj=partstruct{ii}
        partcounter=partcounter+1;
        switchcounter=0;
        for kk=arcstruct{jj}
            ctrp=ctrp+1;
            switchcounter=switchcounter+1;
            ctr=ctr+1;
            if switchcounter>1
                switchtimeindex{ii}=[switchtimeindex{ii} ctr];
            end
        end
        if partcounter>1
            if optflag
                optimalparttimeindex{ii}=[optimalparttimeindex{ii} ctr];
            else
                parttimeindex{ii}=[parttimeindex{ii} ctr];
            end
            switchparttimeindex{ii}=[switchparttimeindex{ii} ctr];
        end
        optflag=OCMATINDIF.optimalswitchingtime(ctrp);
    end
    ctr=ctr+1;
end
OCMATINDIF.switchparttimeindex=switchparttimeindex;
OCMATINDIF.switchtimeindex=switchtimeindex;
OCMATINDIF.parttimeindex=parttimeindex;
OCMATINDIF.optimalparttimeindex=optimalparttimeindex;
if ~isempty(OCMATINDIF.freestatevector)
    OCMATINDIF.freestatevectorcoordinate=length(sol.parameters)+(1:size(OCMATINDIF.freestatevector,2));
    sol.parameters=[sol.parameters zeros(1,size(OCMATINDIF.freestatevector,2))];
end
if ~isempty(freeparametervector)
    OCMATINDIF.freeparametercoordinate=length(sol.parameters)+(1:length(OCMATINDIF.freeparametervectorcoordinate));
    sol.parameters=[sol.parameters zeros(1,length(OCMATINDIF.freeparametervectorcoordinate))];
end
for ii=1:indifforder
    OCMATINDIF.switchtimecoordinate{ii}=length(sol.parameters)+(1:length(switchtimeindex{ii}));
    sol.parameters=[sol.parameters sol.arcinterval(switchtimeindex{ii})];

    OCMATINDIF.parttimecoordinate{ii}=length(sol.parameters)+(1:length(parttimeindex{ii}));
    sol.parameters=[sol.parameters sol.arcinterval(parttimeindex{ii})];

    OCMATINDIF.opttimecoordinate{ii}=length(sol.parameters)+(1:length(optimalparttimeindex{ii}));
    sol.parameters=[sol.parameters sol.arcinterval(optimalparttimeindex{ii})];
end
OCMATINDIF.continuationcoordinate=length(sol.parameters)+1;
sol.parameters=[sol.parameters continuationparameter];


ctr=0;
arcctr=0;
arcctr0=0;
for ii=1:indifforder
    for jj=partstruct{ii}
        ctr=ctr+1;
        OCMATINDIF.edge{ctr}=[OCMATINDIF.arcargument{ii}(1:end-1);OCMATINDIF.arcargument{ii}(2:end)];
        OCMATINDIF.numarc(ctr)=length(OCMATINDIF.arcargument{ctr});
        for kk=1:length(arcstruct{ctr})
            arcctr=arcctr+1;
            OCMATINDIF.solutionindex(arcctr)=ii;
            OCMATINDIF.partindex(arcctr)=jj;
        end
    end
    OCMATINDIF.numarcpersolution(ii)=arcctr-arcctr0+1;
    arcctr0=arcctr;
end
OCMATINDIF.cumsumnumarc=cumsum(OCMATINDIF.numarc);
OCMATINDIF.numarcpersolution=[0 cumsum(OCMATINDIF.numarcpersolution)];

OCMATINDIF.initcoord=[1 OCMATINDIF.cumsumnumarc(1:end-1)+1];
OCMATINDIF.endcoord=OCMATINDIF.cumsumnumarc(1:end);

OCMATINDIF.statecoordinate=scoord;

pathname=OCMATINDIF.datapath{1}();
[resultfile,globalvarfile]=OCMATINDIF.saveintermediatefiles{1}();
OCMATINDIF.basicresultfilename=fullocmatfile(pathname,resultfile);
OCMATINDIF.basicglobalvarfilename=fullocmatfile(pathname,globalvarfile);

OCMATCONT.HE.numinitialcondition=[];

OCMATCONT.codimension=1;




%%%%%%%%%%%%%%%%%%%%%%%%%%

% part denotes in a multistage model a specific stage, the indiffence
% solution is combined to one solution with combined part structure
parts2indifference=zeros(1,numparts); % number of indifference solution for specific part of the path 
parts2numarc=zeros(1,numparts); % number of arcs for a spicific part
arc2part=zeros(1,numparts); % number of arcs for a spicific part
arcindex4optimaltime=zeros(1,numparts); % number of arcs for a spicific part
arc2optimaltimeindex=zeros(1,numparts); % number of arcs for a spicific part
arc2indifference=zeros(1,numparts); % number of arcs for a spicific part
arc2relarc=zeros(1,numparts); % number of arcs for a spicific part
arc2reltimeshift=zeros(1,numparts); % number of arcs for a spicific part
parts2switchtimeindex=cell(1,numparts);
parts2partimeindex=cell(1,numparts);
parts2arccoordinate=cell(1,numparts);
parts2localpartsindex=zeros(1,numparts);
indiffence2arcindex=cell(1,indifforder);
indiffence2partcoordinate=cell(1,indifforder); % the location of the boundary states for the parts among the total number of boundary states for the specific indiffernce solution
indiffence2parttimeindex=cell(1,indifforder); % 
indiffence2partcounter=cell(1,indifforder); % the counter of arcs for each indifference solution, e.g. if in indiff=2, 5 parts exist it is 1:5
indiffence2arcarguments=cell(1,indifforder);

opttimecoordinate=[OCMATINDIF.opttimecoordinate{:}];
%parts2timeindex=cell(1,numparts);
numarccounter=0;
arccounter=0;
timeshift=0;
for ii=1:indifforder
    timecounter=0;
    numarccounter0=numarccounter+1;
    indiffence2partcoordinate{ii}=[];
    numarccounter2=1;
    counter=0;
    relarccounter=0;
    for jj=partstruct{ii}
        counter=counter+1;
        numarccounter1=numarccounter+1;
        timecounter0=timecounter+1;
        arcarg=arcargument(ocTrjMP(jj));
        parts2indifference(jj)=ii;
        parts2numarc(jj)=length(arcarg);
        indiffence2arcarguments{ii}{counter}=arcarg;
        timecounter=timecounter+length(arcarg)+1;
        parts2switchtimeindex{jj}=timecounter0+1:timecounter-1;
        parts2partimeindex{jj}=[timecounter0 timecounter];
        timecounter=timecounter-1;
        numarccounter=numarccounter+length(arcarg);
        parts2arccoordinate{jj}=numarccounter1:numarccounter;
        indiffence2partcoordinate{ii}=[indiffence2partcoordinate{ii} numarccounter];
        numarccounter2=numarccounter2+length(arcarg);
        parts2localpartsindex(jj)=counter;
        flag4opttime=1;
        for kk=1:length(arcarg)
            arccounter=arccounter+1;
            relarccounter=relarccounter+1;
            arc2part(arccounter)=jj;
            arc2indifference(arccounter)=ii;
            arc2relarc(arccounter)=relarccounter;
            arc2reltimeshift(arccounter)=timeshift;
            if OCMATINDIF.optimalswitchingtime(jj) && flag4opttime
                arcindex4optimaltime(arccounter)=1;
                arc2optimaltimeindex(arccounter)=opttimecoordinate(1);
                opttimecoordinate(1)=[];
                flag4opttime=0;
            else
                arcindex4optimaltime(arccounter)=0;
            end
        end
    end
    timeshift=timeshift+numarccounter;
    indiffence2arcindex{ii}=numarccounter0:numarccounter;
    indiffence2parttimeindex{ii}=1:numarccounter2;
    indiffence2partcounter{ii}=1:length(partstruct{ii});
end
OCMATINDIF.parts2indifference=parts2indifference;
OCMATINDIF.parts2numarc=parts2numarc;
OCMATINDIF.parts2switchtimeindex=parts2switchtimeindex;
OCMATINDIF.parts2partimeindex=parts2partimeindex;
OCMATINDIF.parts2arccoordinate=parts2arccoordinate;
OCMATINDIF.indiffence2arcindex=indiffence2arcindex;
OCMATINDIF.indiffence2partcoordinate=indiffence2partcoordinate;
OCMATINDIF.indiffence2parttimeindex=indiffence2parttimeindex;
OCMATINDIF.indiffence2partcounter=indiffence2partcounter;
OCMATINDIF.parts2localpartsindex=parts2localpartsindex;
OCMATINDIF.indiffence2arcarguments=indiffence2arcarguments;
OCMATINDIF.arc2part=arc2part;
OCMATINDIF.arc2indifference=arc2indifference;
OCMATINDIF.arc2relarc=arc2relarc;
OCMATINDIF.arc2reltimeshift=arc2reltimeshift;
OCMATINDIF.arcindex4optimaltime=arcindex4optimaltime;
OCMATINDIF.arc2optimaltimeindex=arc2optimaltimeindex;
if OCMATINDIF.continuationtype==1

    counter=1;
    arccounter=0;
    conttimecounter=0;
    conttimeindex2arc=zeros(1,length(OCMATINDIF.continuationindex));
    for ii=1:indifforder
        for jj=partstruct{ii}
            for kk=1:length(arcarg)
                if any(counter==OCMATINDIF.continuationindex)
                    conttimecounter=conttimecounter+1;
                    conttimeindex2arc(conttimecounter)=arccounter;
                end
                arccounter=arccounter+1;
                counter=counter+1;
            end
        end
    end
    counter=counter+1;
    if any(counter==OCMATINDIF.continuationindex)
        conttimecounter=conttimecounter+1;
        conttimeindex2arc(conttimecounter)=arccounter;
    end
    OCMATINDIF.conttimeindex2arc=conttimeindex2arc;
end
function sol=generatesolstruct(ocMultiPath)
partstruct=partstructure(ocMultiPath);

sol.x=[];
sol.y=[];
sol.arcarg=[];
sol.arcinterval=[];
ctr=0;
for ii=1:length(partstruct)
    for jj=partstruct{ii}
        ctr=ctr+1;
        sol.x=[sol.x independentvar(ocMultiPath(ctr))+ctr-1];
        sol.y=[sol.y dependentvar(ocMultiPath(ctr))];
        sol.arcarg=[sol.arcarg arcargument(ocMultiPath(ctr))];
        sol.arcinterval=[sol.arcinterval arcinterval(ocMultiPath(ctr))];
    end
end
sol.arcinterval(diff(sol.arcinterval)==0)=[];
sol.parameters=[];
sol.x0=sol.arcinterval(1);
sol.solver='';
sol.parameters=[];
sol.solverinfo.tangent=[];
sol.solverinfo.coeff=[];
sol.solverinfo.tmesh=[];
