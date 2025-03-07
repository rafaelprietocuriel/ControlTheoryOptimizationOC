function sol=initocmat_AE_IS4implicit(ocObj,ocMP,varargin)
%
% INITOCMAT_AE_IS initialization for the continuation of an indifference
% threshold
%
% SOL=INITOCMAT_AE_IS(OCOBJ,OCMP,TARGETCOORDINATE,TARGETVALUE) the
% continuation of an indifference threshold in the state space is
% initialized.
% OCOBJ          ... corresponding optimal control model
% OCMP           ... two cell array of ocasymptotics, for the different
%                    solution paths or an instance of an ocmultipath
%                    object.
% TARGETCOORDINATE ... the continuation is done along the n-1 coordinates
%                   (n number of states)
% TARGETVALUE    ... The value of the target vector for the n-1
%                   coordinates.
%
% The output argument SOL is a solution structure for a BVP with the
% standard fields 'x', 'y', 'parameters' (MATLAB syntax for BVPs).
% The (important) OCMat specific fields
%   arcinterval ... truncation of the time interval
%   arcarg      ... arc identifier for a specific combination of active and
%                   inactive constraints.
% are provided as well.
%
% During the initialization two global variables OCMATCONT and OCMATINDIF
% are initialized. The first variable contains general information for the
% continuation process. The second variable contains problem specific
% information.
%
% SOL=INITOCMAT_AE_IS(OCOBJ,OCMP,TARGETCOORDINATE,TARGETVALUE,OPT)

global OCMATCONT OCMATINDIF
sol=[];

% input argument ocMP is either a cell of ocasymptotics or a multi path object
ocMP=ocmultipath(ocMP);
indifforder=multiplicity(ocMP);
targetvalue=[];
stopcriterion=[];
freeendtime=[];
fixdistance=[];
fixinitstate=[];
fixdistancecoordinate=[];
fixendcoordinate=[];
targetcoordinate=[];
freevector=[];
targetvectorcoordinate=[];
asymptoticmatrix=[];
option=[];
divergingcoordinate=[];
exogenousfunction=[];
exogenousinitialstates=[];
userfunction=[];
hitfunction=[];
freeparameter=[];
targetparametervalue=[];
targetparameter=[];
simple=[];
player=[];

for ii=1:indifforder
    pathtpe{ii}=pathtype(ocMP(ii));
end
if isempty(ocObj)
    ocmatmsg('oc model is empty.')
    return
end
if isempty(ocMP)
    ocmatmsg('oc trajectory is empty.')
    return
end
if nargin==2
    option=defaultocoptions;
end

for ii=1:2:length(varargin)
    if ischar(varargin{ii})
        eval([varargin{ii} '=varargin{ii+1};'])
    end
end

if isempty(stopcriterion)
    stopcriterion=1; % stop if limitpoint occurs
end
if isempty(freeendtime)
    freeendtime=zeros(1,indifforder);
end
if isempty(asymptoticmatrix)
    asymptoticmatrix=cell(1,indifforder);
end
OCMATINDIF.stopcriterion=stopcriterion;
if numel(freeendtime)==1
    freeendtime=repmat(freeendtime,1,indifforder);
end
for ii=1:indifforder
    limSet{ii}=limitset(ocMP(ii));
%     if isperiodic(limSet{ii}) && strcmp(pathtpe{ii},'s')
%         freeendtime(ii)=-1;
%     end
end
if isempty(option)
    option=defaultocoptions;
end
if isempty(fixdistance)
    fixdistance=zeros(1,indifforder);
end
if isempty(fixdistancecoordinate)
    fixdistancecoordinate=cell(1,indifforder);
end
if ~isempty(divergingcoordinate)
    if ~iscell(divergingcoordinate)
        divergingcoordinate=repmat({divergingcoordinate},1,indifforder);
    end
end
if ~isempty(fixendcoordinate)
    if ~iscell(fixendcoordinate)
        fixendcoordinate=repmat({fixendcoordinate},1,indifforder);
    end
else
    fixendcoordinate=cell(1,indifforder);
end
if ~isempty(freeparameter)
    freeparameter=parameterindex(ocObj,freeparameter);
end

if isempty(simple)
    simple=zeros(1,indifforder);
end

OCMATINDIF.userfunction=~isempty(userfunction);
OCMATINDIF.hitfunction=~isempty(hitfunction);

OCMATINDIF.exogenousfunction=exogenousfunction;
if OCMATINDIF.exogenousfunction
    OCMATINDIF.exogenousinitialstates=exogenousinitialstates;
    OCMATINDIF.exogenousnumberofstates=length(exogenousinitialstates);
end

if isempty(simple)
    simple=zeros(1,indifforder);
end

ZeroDeviationTolerance=zeros(1,indifforder);
AsymptoticBCMethod=cell(1,indifforder);
for ii=1:indifforder
    if iscell(option)
        ZeroDeviationTolerance(ii)=getocoptions(option{ii},'GENERAL','ZeroDeviationTolerance'); % numeric tolerance for a value to be accepted as zero
        AsymptoticBCMethod{ii}=getocoptions(option{ii},'OCCONTARG','AsymptoticBCMethod');
    else
        ZeroDeviationTolerance(ii)=getocoptions(option,'GENERAL','ZeroDeviationTolerance'); % numeric tolerance for a value to be accepted as zero
        AsymptoticBCMethod{ii}=getocoptions(option,'OCCONTARG','AsymptoticBCMethod');
    end
end
OCMATCONT.ZeroDeviationTolerance=ZeroDeviationTolerance;
OCMATCONT.AsymptoticBCMethod=AsymptoticBCMethod;
% for ii=1:indifforder
%     % remove solverinfo since actual BVP solver and solver used for the
%     % calculation of ocAsym are not identical
%     ocMP(ii).solver='';
%     ocMP(ii).solverinfo=[];
% end
% initialize global variable (OCMATCONT) for general continuation process
OCMATCONT.modelname=modelname(ocObj);
OCMATCONT.modelfunc=modelspecificfunc(ocObj,'4IndifferenceSolutionContinuation');

% initialize global variable (OCMATINDIF) for specific continuation process
% initialize functions
% empty function handles mean either that they have to be computed
% numerically (derivatvies, jacobian, hessian) or that standard functions
% are
% used (e.g. asymptotic boundary condition, initial condition)
funch=OCMATCONT.modelfunc(); % model specific function handles for saddle path continuation

OCMATINDIF.canonicalsystem=funch{1};
OCMATINDIF.canonicalsystemjacobian=funch{2}{1};
OCMATINDIF.canonicalsystemparameterjacobian=funch{2}{2};
OCMATINDIF.canonicalsystemhessian=funch{3}{1};
OCMATINDIF.canonicalsystemparameterhessian=funch{3}{2};
OCMATINDIF.dimplicitcontroldx=funch{4};

% if OCMATINDIF.exogenousfunction
%     OCMATINDIF.exogenousdynamics=funch{4}{1};
%     OCMATINDIF.exogenousjacobian=funch{4}{2};
%     OCMATINDIF.exogenousparameterjacobian=funch{4}{3};
%     try
%         OCMATINDIF.exogenousderivativetime=funch{4}{4};
%     catch
%         OCMATINDIF.exogenousderivativetime=[];
%     end
%     OCMATINDIF.exogenousinitialstatesfunc=funch{4}{11};
% end
% 
% function for the boundary conditions
OCMATINDIF.bcinitial=funch{5}{1};
OCMATINDIF.bcasymptotic=funch{5}{2};
OCMATINDIF.bctransversality=funch{5}{3};
if ~isempty(divergingcoordinate)
    OCMATINDIF.bcinf=funch{5}{4};
else
    OCMATINDIF.equilibrium=funch{5}{4};
end
OCMATINDIF.bcindifference=funch{5}{5};
OCMATINDIF.algebraicequation=funch{5}{7};

% function for Jacobian
OCMATINDIF.bcjacobianinitial=funch{6}{1};
OCMATINDIF.bcjacobianasymptotic=funch{6}{2};
OCMATINDIF.bcjacobiantransversality=funch{6}{3};
OCMATINDIF.bcjacobianindifference=funch{6}{4};

% function describing the hybrid structure of the problem
OCMATINDIF.hybridinfo=funch{7}{1};
OCMATINDIF.domain=funch{7}{2};
OCMATINDIF.guard=funch{7}{3};
OCMATINDIF.reset=funch{7}{4};
%OCMATINDIF.switchtime=funch{7}{5};
%OCMATINDIF.jacobianguard=funch{7}{7};
OCMATINDIF.jacobianreset=funch{7}{8};
OCMATINDIF.domaindiscretization=funch{7}{9};
OCMATINDIF.timesettransformation=funch{7}{10};

if OCMATINDIF.hitfunction
    OCMATINDIF.targetfunction=funch{10}{1};
end
if OCMATINDIF.userfunction
    OCMATINDIF.userfunction=funch{10}{2};
end

% general function
OCMATINDIF.plotcontinuation=funch{11};
OCMATINDIF.testadmissibility=funch{12};
OCMATINDIF.datapath=funch{20};
OCMATINDIF.saveintermediatefiles=funch{21};

hybridinfo=OCMATINDIF.hybridinfo();
% retrieve information about the arcs from the specifc model file
% domaindata and discretizationdata are ordered in the sequence the
% arcarguments appear in hybridinfo.arcarg
for ii=1:numel(hybridinfo.arcarg)
    domaindata(ii)=OCMATINDIF.domain(hybridinfo.arcarg(ii));
end
for ii=1:numel(domaindata)
    OCMATCONT.DOMAINDDATA(ii).numode=domaindata(ii).odedim;
    OCMATCONT.DOMAINDDATA(ii).numae=domaindata(ii).aedim;
    OCMATCONT.DOMAINDDATA(ii).daeorder=domaindata(ii).daeorder;
    OCMATCONT.DOMAINDDATA(ii).numeq=numel(domaindata(ii).daeorder);%number of equations
    OCMATCONT.DOMAINDDATA(ii).eqcoord=1:OCMATCONT.DOMAINDDATA(ii).numeq;
    OCMATCONT.DOMAINDDATA(ii).odecoord=1:domaindata(ii).odedim;
    OCMATCONT.DOMAINDDATA(ii).aecoord=domaindata(ii).odedim+(1:domaindata(ii).aedim);
end
limsetcounter=1;
octrajectory2limset=zeros(indifforder,2);
limSet=cell(1,indifforder);
limitsettype=cell(1,indifforder);
limsetarcarg=cell(1,indifforder);
limsetarcint=cell(1,indifforder);
numberofodes=[];
for ii=1:indifforder
    OCMATINDIF.numarc(ii)=arcnum(ocMP(ii));
    numberofodes=[numberofodes odenumber(ocMP(ii))];

    pathtpe{ii}=pathtype(ocMP(ii));
    OCMATINDIF.stableflag{ii}=~isempty(strfind(pathtpe{ii},'s'));
    if ii>1
        isident=zeros(1,limsetcounter);
        for jj=1:limsetcounter
            isident(jj)=isidentical(limSet{jj},limitset(ocMP(ii)));
        end
        if 1%~all(isident)
            limsetcounter=limsetcounter+1;
            limSet{limsetcounter}=limitset(ocMP(ii));
        end
        octrajectory2limset(ii,:)=[ii limsetcounter];
    else
        octrajectory2limset(1,:)=1;
        limSet{ii}=limitset(ocMP(ii));
    end
    if isequilibrium(limSet{limsetcounter})
        limitsettype{limsetcounter}='e';
    else
        limitsettype{limsetcounter}='l';
    end
    limsetarcarg{limsetcounter}=arcargument(limSet{limsetcounter});
    limsetarcint{limsetcounter}=arcinterval(limSet{limsetcounter});
    OCMATINDIF.implicitcontrolnum(ii)=length(implicitcontrolcoordinate(ocObj,limsetarcarg{limsetcounter}(end)));
end

% model information
OCMATINDIF.statecoordinate=statecoord(ocObj);
OCMATINDIF.statecostatecoordinate=1:2*statenum(ocObj);
OCMATINDIF.player=player;

% test if pure state constraints are defined
OCMATINDIF.stateconstraint=stateconstraint(ocObj);
OCMATINDIF.freeendtime=freeendtime;
OCMATINDIF.fixdistance=fixdistance;
for ii=1:indifforder
    if fixdistance(ii) && isempty(fixdistancecoordinate{ii})
        OCMATINDIF.fixdistancecoordinate{ii}=OCMATINDIF.statecostatecoordinate;
    else
        OCMATINDIF.fixdistancecoordinate{ii}=fixdistancecoordinate{ii};
    end
end
OCMATINDIF.divergingcoordinate=divergingcoordinate;
OCMATINDIF.fixendcoordinate=fixendcoordinate;
OCMATINDIF.fixinitstate=fixinitstate;
OCMATINDIF.freeparameter=freeparameter;
OCMATINDIF.simple=simple;

OCMATINDIF.targetvalue=targetvalue;
OCMATINDIF.targetcoordinate=targetcoordinate;
OCMATINDIF.targetvectorcoordinate=targetvectorcoordinate;

OCMATINDIF.targetparameter=parameterindex(ocObj,targetparameter);
OCMATINDIF.targetparametervalue=targetparametervalue;

OCMATINDIF.indifferenceorder=indifforder;

OCMATINDIF.exogenousfunction=exogenousfunction;
if exogenousfunction
    OCMATINDIF.exogenousinitialstates=exogenousinitialstates;
    OCMATINDIF.exogenousnumberofstates=length(exogenousinitialstates);
end


if ~isempty(divergingcoordinate)
    for ii=1:OCMATINDIF.indifferenceorder
        OCMATINDIF.convergingcoordinate{ii}=setdiff(OCMATINDIF.statecostatecoordinate,OCMATINDIF.divergingcoordinate{ii});
    end
else
     OCMATINDIF.divergingcoordinate=cell(OCMATINDIF.indifferenceorder,1);
end

sol=generatesolstruct(ocMP,limSet);
parametercoordinate=length(sol.parameters);
ctre=0;
ctrl=0;
for ii=1:limsetcounter
    if isequilibrium(limSet{ii})
        ctre=ctre+1;
        % tmp=dependentvar(limSet{ii});
        % limSetdepvar=tmp{1};
        % sol.parameters=[sol.parameters limSetdepvar(:).'];
        % OCMATINDIF.equilibriumcoordinate{ctre}=parametercoordinate+(1:length(limSetdepvar));
        % parametercoordinate=length(sol.parameters);
        octrajectory2limset(octrajectory2limset(:,2)==ii,3)=ctre;
    else
        ctrl=ctrl+1;
        OCMATCONT.monodromy=linearization(limSet{ii});
        dxdt=canonicalsystem(ocObj,limSet{ii},[],1);
        octrajectory2limset(octrajectory2limset(:,2)==ii,3)=ctrl;

        % for the phase condition
        OCMATINDIF.velocityvector{ctrl}=dxdt(:,1);
        OCMATINDIF.velocitycoordinate{ctrl}=1:length(dxdt(:,1));
        OCMATINDIF.velocityvector{ctrl}=OCMATINDIF.velocityvector{ctrl}/norm(OCMATINDIF.velocityvector{ctrl});
        OCMATINDIF.initialpoint{ctrl}=dependentvar(limSet{ii});
        OCMATINDIF.initialpoint{ctrl}=OCMATINDIF.initialpoint{ctrl}(:,1);
    end
end
OCMATINDIF.limsetcounter=limsetcounter;
OCMATINDIF.equilibriumcounter=ctre;
OCMATINDIF.limitcyclecounter=ctrl;
OCMATINDIF.limitsettype=limitsettype;
OCMATINDIF.limsetarcarg=limsetarcarg;
OCMATINDIF.limsetarcint=limsetarcint;
OCMATINDIF.octrajectory2limset=octrajectory2limset;

OCMATINDIF.distance=zeros(1,indifforder);
for ii=1:length(OCMATINDIF.fixdistance)
    tmp=dependentvar(limSet{octrajectory2limset(ii,2)});
    OCMATINDIF.saddlepoint{ii}=tmp{1}(OCMATINDIF.statecostatecoordinate);
    if OCMATINDIF.fixdistance(ii)
        tmp=dependentvar(ocMP(ii));
        OCMATINDIF.distance(ii)=norm(OCMATINDIF.saddlepoint{ii}-tmp{end}(OCMATINDIF.statecostatecoordinate,end));
    end
end
% mode and path specific variables

for ii=1:OCMATINDIF.indifferenceorder
    if isequilibrium(limSet{octrajectory2limset(ii,2)})
        maptype='c';
    else
        maptype='d';
    end
    OCMATINDIF.maptype{ii}=maptype;

    tmp=jacobian(limSet{octrajectory2limset(ii,2)});
    J=tmp{1};
    switch OCMATINDIF.limitsettype{ii}
        case 'e'
            if isempty(asymptoticmatrix{ii})
                Jred=J;
                if ~isempty( OCMATINDIF.divergingcoordinate{ii})
                    Jred(OCMATINDIF.divergingcoordinate{ii},:)=[];
                    Jred(:,OCMATINDIF.divergingcoordinate{ii})=[];
                end
                %[OCMATINDIF.asymptoticmatrix{ii},numstable,numunstable,numcenter,infoStruct]=asymptoticbc(Jred,pathtpe{ii},maptype,ZeroDeviationTolerance(ii));
                [OCMATINDIF.asymptoticmatrix{ii},numstable,numunstable,numcenter,infoStruct]=asymptoticbc(Jred,pathtpe{ii},maptype,ZeroDeviationTolerance(ii),AsymptoticBCMethod{ii});
                if ~isempty(OCMATINDIF.freeparameter) && ~OCMATINDIF.simple(ii)
                    OCMATINDIF.numstable{ii}=numstable;
                    OCMATINDIF.numunstable{ii}=numunstable;
                    OCMATINDIF.numcenter{ii}=numcenter;
                end
            else
                OCMATINDIF.asymptoticmatrix{ii}=asymptoticmatrix{ii};
            end
            tmp=dependentvar(limSet{ii});
            limSetdepvar=tmp{1};
            OCMATINDIF.saddlepoint{ii}=limSetdepvar;
            if ~isempty(OCMATINDIF.freeparameter)
                sol.parameters=[sol.parameters limSetdepvar(:).'];
                OCMATINDIF.equilibriumcoordinate{ii}=parametercoordinate+(1:length(limSetdepvar));
                parametercoordinate=length(sol.parameters);
                if ~OCMATINDIF.simple(ii)
                    switch pathtpe{ii}
                        case 's'
                            OCMATINDIF.subspacedim{ii}=OCMATINDIF.numstable{ii};
                        case {'u','stu'}
                            OCMATINDIF.subspacedim{ii}=OCMATINDIF.numunstable{ii};
                        case {'sc','cs'}
                            OCMATINDIF.subspacedim{ii}=OCMATINDIF.numstable{ii}+OCMATINDIF.numcenter{ii};
                    end
                    OCMATINDIF.orthspacedim{ii}=OCMATINDIF.statecostatecoordinate(end)-OCMATINDIF.subspacedim{ii};
                    Y=zeros(OCMATINDIF.orthspacedim{ii},OCMATINDIF.subspacedim{ii});
                    OCMATINDIF.Y{ii}=Y;
                    OCMATINDIF.Q0{ii}=infoStruct.Q;
                    OCMATINDIF.Ycoordinate{ii}=reshape(parametercoordinate+(1:OCMATINDIF.orthspacedim{ii}*OCMATINDIF.subspacedim{ii}),OCMATINDIF.orthspacedim{ii},OCMATINDIF.subspacedim{ii});
                    sol.parameters=[sol.parameters Y(:).'];
                    parametercoordinate=length(sol.parameters);
                    OCMATINDIF.Id{ii}=eye(OCMATINDIF.orthspacedim{ii});
                    OCMATINDIF.numY{ii}=numel(Y);
                end
            end
            if ~isempty(OCMATINDIF.fixendcoordinate{ii})
                depvar=dependentvar(ocMP(ii));
                OCMATINDIF.endvalue{ii}=depvar(OCMATINDIF.fixendcoordinate{ii},end);
            else
                OCMATINDIF.endvalue{ii}=[];
            end
        case 'l'
    end
end
depvar=dependentvar(ocMP(1));

if OCMATINDIF.exogenousfunction
    OCMATINDIF.exogenousdynamicscoordinate=numberofodes+1:numberofodes+OCMATINDIF.exogenousnumberofstates;
else
    OCMATINDIF.exogenousdynamicscoordinate=[];
end

OCMATINDIF.cumsumnumarc=cumsum(OCMATINDIF.numarc);
OCMATINDIF.initcoord=[1 OCMATINDIF.cumsumnumarc(1:end-1)+1];

OCMATINDIF.parametervalue=parametervalue(ocObj);
OCMATINDIF.initialtime=sol.x0;
%lefttimeindex=[1 cumsum(OCMATINDIF.numarc(1:end-1)+1)];
if ~isempty(OCMATINDIF.freeparameter)
    OCMATINDIF.freeparametercoordinate=parametercoordinate+(1:length(OCMATINDIF.freeparameter));
    sol.parameters=[sol.parameters OCMATINDIF.parametervalue(OCMATINDIF.freeparameter)];
    parametercoordinate=length(sol.parameters);
end
OCMATINDIF.freevectorcoordinate=parametercoordinate+(1:size(freevector,2));
for ii=1:size(freevector,2)
    sol.parameters=[sol.parameters 0];
end
OCMATINDIF.freevector=freevector;
OCMATINDIF.linearization=J;
OCMATINDIF.limitset=limSet;
OCMATINDIF.startvalue=depvar{1}(OCMATINDIF.statecoordinate,1);
OCMATINDIF.pathtype=pathtpe;
pathname=OCMATINDIF.datapath();
[resultfile,globalvarfile]=OCMATINDIF.saveintermediatefiles();
OCMATINDIF.basicresultfilename=fullocmatfile(pathname,resultfile);
OCMATINDIF.basicglobalvarfilename=fullocmatfile(pathname,globalvarfile);
OCMATINDIF.objectivevaluecalc=0;
OCMATINDIF.autonomous=1;
OCMATINDIF.objectivevaluecoord=[];


OCMATCONT.HE.numinitialcondition=numel(targetvalue);
OCMATCONT.HE.numendcondition=size(OCMATINDIF.asymptoticmatrix{1},2)*OCMATINDIF.indifferenceorder;
OCMATCONT.numode=numberofodes;
OCMATCONT.maxnumode=max(numberofodes);

OCMATCONT.codimension=1;

function sol=generatesolstruct(ocMultiPath,limSet)
global OCMATINDIF
OCMATINDIF.solutionindex=[];
fn={'x','y','arcarg','arcinterval'};
ctr=0;
parametercoordinate=0;
for ii=1:OCMATINDIF.indifferenceorder
    if ii==1
        sol=odestruct(ocMultiPath(1));
        maxodenum=size(sol.y,1);
    else
        soltmp=odestruct(ocMultiPath(ii));
        maxodenum=max([maxodenum size(soltmp.y,1)]);
        y=sol.y;
        sol.y=zeros(maxodenum,size(y,2));
        sol.y(1:size(y,1),:)=y;
        y=soltmp.y;
        soltmp.y=zeros(maxodenum,size(y,2));
        soltmp.y(1:size(y,1),:)=y;
        soltmp.x=soltmp.x+sol.x(end);
        soltmp.arcposition=soltmp.arcposition+length(sol.x);
        for jj=1:length(fn)
            sol.(fn{jj})=[sol.(fn{jj}) soltmp.(fn{jj})];
        end
        sol.arcpositon=[sol.arcposition soltmp.arcposition];
    end
    ctr=ctr+1;
    arcintv=arcinterval(ocMultiPath(ii));
    arcarg=arcargument(ocMultiPath(ii));


    OCMATINDIF.switchtime{ctr}=arcintv(2:end-1);
    OCMATINDIF.truncationtime(ctr)=arcintv(end);
    OCMATINDIF.numarc(ctr)=length(arcarg);
    OCMATINDIF.arcarg{ctr}=arcarg;
    OCMATINDIF.edge{ctr}=[arcarg(1:end-1);arcarg(2:end)];
    OCMATINDIF.solutionindex=[OCMATINDIF.solutionindex repmat(ctr,1,OCMATINDIF.numarc(ctr))];

    sol.parameters=[sol.parameters OCMATINDIF.switchtime{ctr}];
    OCMATINDIF.switchtimecoordinate{ctr}=parametercoordinate+(1:length(OCMATINDIF.switchtime{ctr}));
    parametercoordinate=length(sol.parameters);
    if OCMATINDIF.freeendtime(ctr)
        sol.parameters=[sol.parameters OCMATINDIF.truncationtime(ctr)];
        OCMATINDIF.truncationtimecoordinate(ctr)=parametercoordinate+1;
        parametercoordinate=length(sol.parameters);
    end
end

ctr2=0;
OCMATINDIF.limitcycleindex=zeros(1,length(OCMATINDIF.solutionindex));
for ii=1:length(limSet)
    if ~isequilibrium(limSet{ii})
        ctr2=ctr2+1;
        ctr=ctr+1;
        arcintv=arcinterval(limSet{ii});
        arcarg=arcargument(limSet{ii});
        sol.x=[sol.x independentvar(limSet{ii})+sol.x(end)];
        sol.y=[sol.y dependentvar(limSet{ii})];
        sol.arcarg=[sol.arcarg arcargument(limSet{ii})];
        sol.arcinterval=[sol.arcinterval arcintv];
        sol.x0(ctr)=0;
        
        OCMATINDIF.switchtime{ctr}=arcintv(2:end-1);
        OCMATINDIF.numarc(ctr)=length(arcarg);
        OCMATINDIF.arcarg{ctr}=arcarg;
        OCMATINDIF.edge{ctr}=[arcarg(1:end-1);arcarg(2:end)];
        OCMATINDIF.solutionindex=[OCMATINDIF.solutionindex repmat(ctr,1,OCMATINDIF.numarc(ctr))];
        OCMATINDIF.limitcycleindex=[OCMATINDIF.limitcycleindex repmat(ctr2,1,OCMATINDIF.numarc(ctr))];

        sol.parameters=[sol.parameters OCMATINDIF.switchtime{ctr}];
        OCMATINDIF.switchtimecoordinate{ctr}=parametercoordinate+(1:length(OCMATINDIF.switchtime{ctr}));
        parametercoordinate=length(sol.parameters);
        
        sol.parameters=[sol.parameters arcintv(end)];
        OCMATINDIF.periodcoordinate(ctr2)=parametercoordinate+1;
        parametercoordinate=length(sol.parameters);
        
    end
end

arcposition=find(diff(sol.x)==0);
sol.arcposition=[[1 arcposition+1];[arcposition numel(sol.x)]];
sol.solver='';
sol.solverinfo.tangent=[];
sol.solverinfo.coeff=[];
sol.solverinfo.tmesh=[];