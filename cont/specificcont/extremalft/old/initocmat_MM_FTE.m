function sol=initocmat_MM_FTE(mmObj,ocTrj,continuationtype,continuationindex,continuationtarget,varargin)
%
% initocmat_FTE_T initialization for asymptotic extremal calculation
%
% SOL=initocmat_FTE_T(OCOBJ,OCASYM,CONTCOORDINATE,TARGETVALUE) the
% continuation process is initialized (started) at a solution OCASYM
% (instance of the ocasymptotic class) that (usually) has been calculated
% during a previous continuation. Therefore, neither the 'TruncationTime' ,
% nor the 'PathType' (see initocmat_FTE_T) can be provided. This
% information is taken from the OCASYM object. 
%
% For a more detailed discussion see INITOCMAT_AE_EP
%
clear global OCMATCONT OCMATFTE
global OCMATCONT OCMATFTE
sol=[];
initialstate=[];
initialcostate=[];
initialendtime=[];
initialarcargument=[];
findoptimalswitchingtime=[];
fixendstate=[];
fixendcostate=[];
fixinitialstate=[];
fixinitcostate=[];
freeparameter=[];
findoptimalparameter=[];

objectivevaluecalc=[];
addpart=[];
opt=[];
if isempty(mmObj)
    ocmatmsg('oc model is empty.')
    return
end

% contiuationtype: time, state, parameter
% continuationindex: (time) number of switching/endtime, (state)
% coordinate, (parameter) index or name of parameter
% continuationtarget: value
optionidx=find(strcmpi(varargin,'option'));
initialstateidx=find(strcmpi(varargin,'initialstate'));
addpartidx=find(strcmpi(varargin,'addpart'));
initialcostateidx=find(strcmpi(varargin,'initialcostate'));
initialendtimeidx=find(strcmpi(varargin,'initialendtime'));
initialarcargumentidx=find(strcmpi(varargin,'initialarcargument'));
fixinitialstateidx=find(strcmpi(varargin,'fixinitialstate'));
fixendstateidx=find(strcmpi(varargin,'fixendstate'));
freeparameteridx=find(strcmpi(varargin,'freeparameter'));
findoptimalswitchingtimeidx=find(strcmpi(varargin,'findoptimalswitchingtime'));
findoptimalparameteridx=find(strcmpi(varargin,'findoptimalparameter'));
objectivevaluecalcidx=find(strcmpi(varargin,'objectivevaluecalc'));
if ~isempty(addpartidx)
    addpart=varargin{addpartidx+1};
end
if ~isempty(fixendstateidx)
    fixendstate=varargin{fixendstateidx+1};
end
if ~isempty(fixinitialstateidx)
    fixinitialstate=varargin{fixinitialstateidx+1};
end
if ~isempty(freeparameteridx)
    freeparameter=varargin{freeparameteridx+1};
end
if ~isempty(findoptimalparameteridx)
    findoptimalparameter=varargin{findoptimalparameteridx+1};
end
if ~isempty(objectivevaluecalcidx)
    objectivevaluecalc=varargin{objectivevaluecalcidx+1};
end
if ~isempty(findoptimalswitchingtimeidx)
    findoptimalswitchingtime=varargin{findoptimalswitchingtimeidx+1};
end
if isempty(objectivevaluecalc)
    objectivevaluecalc=0;
end
if isempty(findoptimalswitchingtime)
    findoptimalswitchingtime=false;
end
if isempty(findoptimalparameter)
    findoptimalparameter=false;
end
OCMATFTE.objectivevaluecalc=objectivevaluecalc;

nummod=numberofmodels(mmObj);
OCMATFTE.numberofmodels=nummod;
if ~iscell(ocTrj)
    OCMATFTE.optimalswitchingtime=zeros(1,nummod);
    if isoctrajectory(ocTrj)
        ocTrjCell{1}=ocTrj;
    elseif ismmultipath(ocTrj)
        OCMATFTE.optimalswitchingtime=optimaltransition(ocTrj);
        for ii=1:numberofparts(ocTrj)
            ocTrjCell{ii}=ocTrj(ii);
        end
    end
else
    ocTrjCell=ocTrj;
    OCMATFTE.optimalswitchingtime=zeros(1,nummod);
end
switch continuationtype
    case {'initialstate','endstate'}
        OCMATFTE.continuationtype=0;

    case {'switchingtime','time'}
        OCMATFTE.continuationtype=1;
    case 'parameter'
        OCMATFTE.continuationtype=2;
        parindex=continuationindex;
        if ischar(parindex)
            parindex=parameterindex(mmObj,parindex);
        end
        OCMATFTE.parameterindex=parindex;
end
if freeparameter
    parindex=varargin{freeparameteridx+1};
    if ischar(parindex)
        parindex=parameterindex(mmObj,parindex);
    end
    OCMATFTE.parameterindex=parindex;
end
OCMATFTE.freeparameter=freeparameter;
OCMATFTE.findoptimalparameter=findoptimalparameter;

%%% generate initial octrajectory if not provided
if addpart
    if (length(ocTrjCell)==nummod-1 || (length(ocTrjCell)==nummod && isempty(ocTrjCell{nummod})))
        if ~isempty(initialstateidx)
            initialstate=varargin{initialstateidx+1};
            initialstate=initialstate(:);
        end
        if ~isempty(initialcostateidx)
            initialcostate=varargin{initialcostateidx+1};
            initialcostate=initialcostate(:);
        end
        if ~isempty(initialendtimeidx)
            initialendtime=varargin{initialendtimeidx+1};
        end
        if ~isempty(initialarcargumentidx)
            initialarcargument=varargin{initialarcargumentidx+1};
        end
        if ~isempty(optionidx)
            opt=varargin{optionidx+1};
        end
        if isempty(opt)
            opt=defaultocoptions;
        end
        if isempty(initialarcargument)
            initialarcargument=0;
        end
        if isempty(initialendtime)
            initialendtime=0;
        end
        if isempty(initialstate)
            initialstate=state(mmObj(nummod-1),ocTrjCell{nummod-1},1);
            initialstate=initialstate(:,end);
        end
        n=getocoptions(opt,'GENERAL','TrivialArcMeshNum');
        ocTrjN.x=linspace(0,1,n);
        initocPt.y=initialstate;
        t=time(mmObj(nummod-1),ocTrjCell{nummod-1},1);
        initocPt.x=t(end);
        initocPt.arcarg=initialarcargument;
        initocPt.arcinterval=[t(end) t(end)];
        initocPt.arcposition=[1;1];
        if isempty(initialcostate)
            initialcostate=costate(mmObj(nummod-1),ocTrjCell{nummod-1},1);
            initialcostate=initialcostate(:,end);
        end
        ocTrjN.y=[initialstate(:);initialcostate];
        ocTrjN.y=ocTrjN.y(:,ones(1,n));
        ocTrjN.arcposition=[1;n];
        ocTrjN.arcinterval=[t(end) t(end)];
        ocTrjN.arcarg=initialarcargument;
        ocTrjN.x0=t(end);
        ocTrjN.timehorizon=t(end);
        ocTrjN.modelparameter=parametervalue(mmObj(nummod));
        ocTrjN.modelname=modelname(mmObj(nummod));
        ocTrjCell{nummod}=octrajectory(ocTrjN);
    end
end

if ~isempty(ocTrjCell)
    ocTrjMP=mmultipath(ocTrjCell);
else
    ocTrjMP=ocTrj;
end

if objectivevaluecalc
    Oval=0;
    for ii=1:nummod
        ocTrjStruct=struct(ocTrjMP(ii));
        if length(ocTrjMP(ii).y(:,1))==2*statenum(mmObj(ii))
            
            o=objectivefunction(mmObj(ii),ocTrjMP(ii),1);
            ocTrjStruct.y(end+1,:)=Oval+[0 cumsum((o(1:end-1)+o(2:end))/2.*diff(time(mmObj(ii),ocTrjMP(ii),1)))];
            ocTrjCell{ii}=octrajectory(ocTrjStruct);
            Oval=ocTrjStruct.y(end,end);
        else
            ocTrjCell{ii}=ocTrjMP(ii);
            Oval=ocTrjMP(ii).y(end,end);

        end
    end
    ocTrjMP=mmultipath(ocTrjCell);
end
scoord=statecoord(mmObj);
cscoord=costatecoord(mmObj);
for ii=1:nummod
    OCMATFTE.statecostatecoord{ii}=[scoord{ii}(:).' cscoord{ii}(:).'];
end

% initialize global variable (OCMATCONT) for general continuation process
OCMATCONT.modelname=modelname(mmObj);
OCMATCONT.modelfunc=modelspecificfunc(mmObj,'4FiniteHorizonPathContinuation');

% initialize global variable (OCMATFTE) for specific continuation process 
% initialize functions 
% empty function handles mean either that they have to be computed
% numerically (derivatives, jacobian, hessian) or that standard functions are
% used (e.g. asymptotic boundary condition, initial condition)
% mode and path specific variables
OCMATFTE.parametervalue=parametervalue(mmObj);
OCMATFTE.autonomous=isautonomous(mmObj);



for ii=1:nummod
    funch{ii}=OCMATCONT.modelfunc{ii}(); % model specific function handles for saddle path continuation
end
for ii=1:nummod
    
    
    OCMATFTE.canonicalsystem{ii}=funch{ii}{1};
    OCMATFTE.canonicalsystemjacobian{ii}=funch{ii}{2}{1};
    OCMATFTE.canonicalsystemparameterjacobian{ii}=funch{ii}{2}{2};
    OCMATFTE.canonicalsystemhessian{ii}=funch{ii}{3}{1};
    OCMATFTE.canonicalsystemparameterhessian{ii}=funch{ii}{3}{2};

    % function for the boundary conditions
    OCMATFTE.bcinitial{ii}=funch{ii}{5}{1};
    OCMATFTE.bctransversality{ii}=funch{ii}{5}{2};
    OCMATFTE.bcoptimalhorizon{ii}=funch{ii}{5}{3};
    OCMATFTE.salvagevalue{ii}=funch{ii}{5}{6};
    OCMATFTE.bcconnectingparts{ii}=funch{ii}{5}{7};
    OCMATFTE.bcoptimalconnectingparts{ii}=funch{ii}{5}{8};
    
    % function for Jacobian
    OCMATFTE.bcjacobianinitial{ii}=funch{ii}{6}{1};
    OCMATFTE.bcjacobiantransversality{ii}=funch{ii}{6}{2};

    % function describing the hybrid structure of the problem
    OCMATFTE.hybridinfo{ii}=funch{ii}{7}{1};
    OCMATFTE.domain{ii}=funch{ii}{7}{2};
    OCMATFTE.guard{ii}=funch{ii}{7}{3};
    OCMATFTE.reset{ii}=funch{ii}{7}{4};
    OCMATFTE.switchtime{ii}=funch{ii}{7}{5};
    OCMATFTE.jacobianguard{ii}=funch{ii}{7}{7};
    OCMATFTE.jacobianreset{ii}=funch{ii}{7}{8};
    OCMATFTE.domaindiscretization{ii}=funch{ii}{7}{9};

    OCMATFTE.objectivefunction{ii}=funch{ii}{8}{1};
    OCMATFTE.objectivefunctionjacobian{ii}=funch{ii}{8}{2};
    OCMATFTE.objectivefunctionparameterjacobian{ii}=funch{ii}{8}{3};
    OCMATFTE.objectivefunctionderivativetime{ii}=funch{ii}{8}{4};

    if ~isautonomous(mmObj(ii))
        OCMATFTE.canonicalsystemderivativetime{ii}=funch{ii}{2}{3};
    end
    % general function
    OCMATFTE.plotcontinuation{ii}=funch{ii}{11};
    OCMATFTE.testadmissibility{ii}=funch{ii}{12};
    OCMATFTE.datapath{ii}=funch{ii}{20};
    OCMATFTE.saveintermediatefiles{ii}=funch{ii}{21};

    hybridinfo{ii}=OCMATFTE.hybridinfo{ii}();
    
    for jj=1:numel(hybridinfo{ii}.arcarg)
        domaindata{ii}(jj)=OCMATFTE.domain{ii}(hybridinfo{ii}.arcarg(jj));
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

OCMATFTE.findoptimalswitchingtime=findoptimalswitchingtime;

sol=generatesolstruct(ocTrjMP);
if  freeparameter
    sol.parameters=[sol.parameters OCMATFTE.parametervalue{1}(parindex{1})];
    OCMATFTE.parametercoord=length(sol.parameters)-length(parindex{1})+(1:length(parindex{1}));
end
depvarFirst=dependentvar(ocTrjMP(1));
depvarLast=dependentvar(ocTrjMP(nummod));
OCMATFTE.initialstate=depvarFirst(scoord{1},1);
OCMATFTE.endstate=depvarLast(scoord{nummod},end);


OCMATCONT.HE.edge=[];
arcoffset=0;
switchtimes=0;
arcintv=arcinterval(ocTrjMP);
arcarg=arcargument(ocTrjMP);
arcn=arcnum(ocTrjMP);
OCMATFTE.partendtime=zeros(1,nummod);
for ii=1:nummod
    switchtimes=[switchtimes arcintv{ii}(2:end)];
    if ii<nummod
        fixparttimeindex(ii)=length(switchtimes);
    end
    OCMATFTE.partendtime(ii)=switchtimes(end);
    OCMATFTE.edge{ii}=[arcarg{ii}(1:end-1);arcarg{ii}(2:end)];
    OCMATFTE.numarc(ii)=arcn{ii};
    OCMATFTE.initialstateindex(ii)=numel(ocTrjMP(ii).x);
    OCMATFTE.arccoord{ii}=(1:arcn{ii})+arcoffset;
    arcoffset=arcoffset+arcn{ii};
end
OCMATFTE.cumsumnumarc=cumsum(OCMATFTE.numarc);
totaltimeindex=1:OCMATFTE.cumsumnumarc(end)+1;
switchtimeindex=setdiff(totaltimeindex,fixparttimeindex); % time indices within one part/stage + initial and end time index
switchtimeindex([1 end])=[];
optimalparttimeindex=fixparttimeindex(OCMATFTE.optimalswitchingtime==1);
counter=0;
for ii=1:nummod
    counter_start=counter+1;
    counter=counter+OCMATFTE.numarc(ii);
    OCMATFTE.solutionindex(counter_start:counter)=ii;
end
OCMATFTE.initialtimeinterval=switchtimes;
switch continuationtype
    case 'initialstate'
        OCMATFTE.continuationvector=continuationtarget(:)-OCMATFTE.initialstate;
        continuationparameter=0;
    case 'endstate'
        OCMATFTE.continuationvector=continuationtarget(:)-OCMATFTE.endstate;
        continuationparameter=0;

    case {'switchingtime','time'}
        OCMATFTE.initialswitchtime=switchtimes(continuationindex);
        OCMATFTE.continuationtimeindex=continuationindex;
        OCMATFTE.continuationvector=continuationtarget(:).'-OCMATFTE.initialswitchtime;
        continuationparameter=0;
    case 'parameter'
        if length(parindex{1})==1
            for ii=1:nummod
                OCMATFTE.initialparameter(ii)=OCMATFTE.parametervalue{ii}(parindex{ii});
            end
            OCMATFTE.continuationvector=continuationtarget(:)-OCMATFTE.initialparameter(:);
            OCMATFTE.initialparameter=OCMATFTE.initialparameter(:);
        else
            OCMATFTE.initialparameter=OCMATFTE.parametervalue{1}(parindex{1});
            OCMATFTE.continuationvector=continuationtarget(:)-OCMATFTE.initialparameter(:);
        end
        continuationparameter=0;
end
OCMATFTE.switchtimeindex=switchtimeindex;
OCMATFTE.switchtimecoordinate=length(sol.parameters)+(1:length(OCMATFTE.switchtimeindex));
sol.parameters=[sol.parameters switchtimes(OCMATFTE.switchtimeindex)];
OCMATFTE.optimalparttimeindex=optimalparttimeindex;
OCMATFTE.optimalpartcoordinate=length(sol.parameters)+(1:length(OCMATFTE.optimalparttimeindex));
sol.parameters=[sol.parameters switchtimes(OCMATFTE.optimalparttimeindex)];
sol.parameters=[sol.parameters continuationparameter];

OCMATFTE.continuationcoordinate=length(sol.parameters);
OCMATFTE.continuationtarget=continuationtarget;
OCMATFTE.continuationindex=continuationindex;


OCMATFTE.initcoord=[1 OCMATFTE.cumsumnumarc(1:end-1)+1];
OCMATFTE.endcoord=OCMATFTE.cumsumnumarc(1:end);


if objectivevaluecalc
    OCMATFTE.objectivevaluecoord=size(sol.y,1);
end
OCMATFTE.statecoordinate=scoord;
OCMATFTE.fixinitialstatecoord=fixinitialstate;
OCMATFTE.fixendstatecoord=fixendstate;
OCMATFTE.fixinitcostatecoord=fixinitcostate;
OCMATFTE.fixendcostatecoord=fixendcostate;
pathname=OCMATFTE.datapath{1}();
[resultfile,globalvarfile]=OCMATFTE.saveintermediatefiles{1}();
OCMATFTE.basicresultfilename=fullocmatfile(pathname,resultfile);
OCMATFTE.basicglobalvarfilename=fullocmatfile(pathname,globalvarfile);

OCMATCONT.HE.numinitialcondition=[];

OCMATCONT.codimension=1;

function sol=generatesolstruct(ocMultiPath)
global OCMATFTE

nummult=numberofparts(ocMultiPath);
sol.x=independentvar(ocMultiPath(1));
sol.y=dependentvar(ocMultiPath(1));
sol.arcarg=arcargument(ocMultiPath(1));
sol.arcinterval=arcinterval(ocMultiPath(1));
sol.parameters=parameters(ocMultiPath(1));
x0=initialtime(ocMultiPath(1));
for ii=2:nummult
    sol.x=[sol.x independentvar(ocMultiPath(ii))-x0+sol.x(end)];

    sol.y=[sol.y dependentvar(ocMultiPath(ii))];
    sol.arcarg=[sol.arcarg arcargument(ocMultiPath(ii))];
    actarcinterval=arcinterval(ocMultiPath(ii));
    sol.arcinterval=[sol.arcinterval actarcinterval(2:end)];
end
sol.parameters=[];
sol.x0=x0;
sol.solver='';
sol.parameters=[];
sol.solverinfo.tangent=[];
sol.solverinfo.coeff=[];
sol.solverinfo.tmesh=[];
