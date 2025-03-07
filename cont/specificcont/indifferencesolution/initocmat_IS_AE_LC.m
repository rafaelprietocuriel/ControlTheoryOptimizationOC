function sol=initocmat_IS_AE_LC(ocObj,ocMP,varargin)
%
% initocmat_AE_PER initialization for asymptotic extremal calculation
%
% SOL=initocmat_AE_PER(OCOBJ,OCEP,CONTCOORDINATE,TARGETVALUE) a stable
% saddle path calculation is initialized.
% OCOBJ          ... corresponding optimal control model
% OCEP           ... equilibrium point (dynprimitive) (hat-x) with (local)
%                    stable manifold of dimension k
% CONTCOORDINATE ... coordinates i_1,...,i_k of the continuation variable
%                    (usually state coordinate(s) in optimal control
%                    problems)
% TARGETVALUE    ... determines direction of the continuation (x_j^0,
%                    j=i_1,...,i_k)
%
% The output argument SOL is a solution structure for a BVP with the
% standard fields 'x', 'y', 'parameters' (MATLAB syntax for BVPs).
% The (important) OCMat specific fields
%   arcinterval ... truncation of the time interval
%   arcarg      ... arc identifier for a specific combination of active and
%                   inactive constraints.
% are provided as well.
%
% During the initialization two global variables OCMATCONT and OCMATINDIF are
% initialized. The first variable contains general information for the
% continuation process. The second variable contains problem specific
% information.
%
% SHORT EXPLANATION OF THE MATHEMATICAL BACKGROUND
%
% The underlying mathematical problem formulation is to find a trajectory
% (time path) x(t)=(x_1(t),...,x_N(t)) with initial condition x_j(0)=x_j^0,
% j=i_1,...,i_k and lim_{t\to\infty}x(t)=hat-x (convergence to the
% equilibrium).
%
% The problem is solved using a continuation algorithm, where the
% continuation is done for the initial condition
%   x_j(0)=x_j^0, j=i_1,...,i_k
% and the continuation parameter 'mu' is defined as
%   x_j(0)=x_j^0*mu+(1-mu)*hat-x_j, j=i_1,...,i_k.
% Thus, for mu=0 we have
%   x_j(0)=hat-x_j, j=i_1,...,i_k
% and for mu=1
%   x_j(0)=x_j^0
% The end condition, convergence to the equilibrium, is reformulated in a
% way that allows a numerical treatment. The default way is the truncation
% of the infinite time to a finite time 'T' and the condtion that the end
% point x(T) ends on the linearized stable manifold (stable eigenspace).
%
% This means that at the start of the continuation the equilibrium path
% (constant solution at the equilibrium) trivially satisfies the initial
% and end condition. Therefore, with the provision of the equilibrium the
% initial solution is given as well. In that sense the OCEP argument
% performs two tasks. The searched for solution converges to OCEP and OCEP
% is the initial solution (mu=0) of the continuation process.
%
% The denomination as TARGETVALUE maybe misleading but from the
% continuation point of view it denotes the target. From the problem
% perspective it denotes the initial state values of the searched for
% solution.
%
% SOL=INITOCMAT_AE_EP(OCOBJ,OCEP,CONTCOORDINATE,TARGETVALUE,OPT) with the
% option structure OPT the threshold 'ZeroDeviationTolerance' and initial
% number of discretization points 'TrivialArcMeshNum' for the equilibrium
% solution can be changed.
%   ZeroDeviationTolerance ... provides the tolerance to classify an
%                              eigenvalue numerically as zero.
%   TrivialArcMeshNum      ... provides the number of points for the
%                              constant solution at the equilbrium.
%
% SOL=INITOCMAT_AE_EP(...,'TruncationTime',T) the truncation of the
% infinite time horizon to the finite time T
%
% SOL=INITOCMAT_AE_EP(...,'PathType',p)
%   p='s' (default) stable saddle-path calculation
%   p='u' unstable saddle-path calculation

clear global OCMATCONT OCMATINDIF
global OCMATCONT OCMATINDIF
sol=[];

ocMP=ocmultipath(ocMP);
indifforder=multiplicity(ocMP);
targetvalue=[];
freeendtime=[];
fixdistance=[];
fixdistancecoordinate=[];
fixinitstate=[];
fixendcoordinate=[];
targetcoordinate=[];
freevector=[];
targetvectorcoordinate=[];
asymptoticmatrix=[];
option=[];
userfunction=[];
hitfunction=[];
freeparameter=[];
targetparametervalue=[];
targetparameter=[];
simple=[];
fixlccoordinate=[];
indifferencevalue=[];
if isempty(ocObj)
    ocmatmsg('oc model is empty.\n')
    return
end

pathtpe=cell(1,indifforder);
for ii=1:indifforder
    if isempty(ocMP(ii))
        return
    end
end


for ii=1:2:length(varargin)
    if ischar(varargin{ii})
        eval([varargin{ii} '=varargin{ii+1};'])
    end
end
if isempty(freeendtime)
    freeendtime=zeros(1,indifforder);
end
if isempty(asymptoticmatrix)
    asymptoticmatrix=cell(1,indifforder);
end
if numel(freeendtime)==1
    freeendtime=repmat(freeendtime,1,indifforder);
end

if isempty(option)
    option=defaultocoptions;
end
if isempty(fixdistance)
    fixdistance=zeros(1,indifforder);
end
if ~isempty(fixendcoordinate)
    if ~iscell(fixendcoordinate)
        fixendcoordinate=repmat({fixendcoordinate},1,indifforder);
    end
else
    fixendcoordinate=cell(1,indifforder);
end
if isempty(fixdistancecoordinate)
    fixdistancecoordinate=cell(1,indifforder);
end
if ~isempty(freeparameter)
    freeparameter=parameterindex(ocObj,freeparameter);
end

if isempty(simple)
    simple=zeros(1,indifforder);
end
if length(simple)==1
    simple=repmat(simple,1,indifforder);
end
if isempty(fixlccoordinate)
    fixlccoordinate=zeros(1,indifforder);
end
if length(fixlccoordinate)==1
    fixlccoordinate=repmat(fixlccoordinate,1,indifforder);
end
OCMATINDIF.fixlccoordinate=fixlccoordinate;

OCMATINDIF.userfunction=~isempty(userfunction);
OCMATINDIF.hitfunction=~isempty(hitfunction);

ZeroDeviationTolerance=zeros(1,indifforder);
asymptoticbcmethod=cell(1,indifforder);
for ii=1:indifforder
    if iscell(option)
        ZeroDeviationTolerance(ii)=getocoptions(option{ii},'GENERAL','ZeroDeviationTolerance'); % numeric tolerance for a value to be accepted as zero
        asymptoticbcmethod{ii}=getocoptions(option{ii},'OCCONTARG','AsymptoticBCMethod');
    else
        ZeroDeviationTolerance(ii)=getocoptions(option,'GENERAL','ZeroDeviationTolerance'); % numeric tolerance for a value to be accepted as zero
        asymptoticbcmethod{ii}=getocoptions(option,'OCCONTARG','AsymptoticBCMethod');
    end
end
% initialize global variable (OCMATCONT) for general continuation process
OCMATCONT.modelname=modelname(ocObj);
OCMATCONT.modelfunc=modelspecificfunc(ocObj,'4IndifferenceSolutionContinuation');

% initialize global variable (OCMATINDIF) for specific continuation process
% initialize functions
% empty function handles mean either that they have to be computed
% numerically (derivatvies, jacobian, hessian) or that standard functions are
% used (e.g. asymptotic boundary condition, initial condition)

funch=OCMATCONT.modelfunc(); % model specific function handles for saddle path continuation

OCMATINDIF.canonicalsystem=funch{1};
OCMATINDIF.canonicalsystemjacobian=funch{2}{1};
OCMATINDIF.canonicalsystemparameterjacobian=funch{2}{2};
OCMATINDIF.canonicalsystemhessian=funch{3}{1};
OCMATINDIF.canonicalsystemparameterhessian=funch{3}{2};


if OCMATINDIF.userfunction
    OCMATINDIF.userfunction=funch{10}{2};
end
if OCMATINDIF.hitfunction
    OCMATINDIF.userfunction=funch{10}{2};
end

% function for the boundary conditions
OCMATINDIF.bcinitial=funch{5}{1};
OCMATINDIF.bcasymptotic=funch{5}{2};
OCMATINDIF.bctransversality=funch{5}{3};
OCMATINDIF.equilibrium=funch{5}{4};
OCMATINDIF.bcindifference=funch{5}{5};
OCMATINDIF.bclimitcycle=funch{5}{6};

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
OCMATINDIF.jacobianguard=funch{7}{7};
OCMATINDIF.jacobianreset=funch{7}{8};
OCMATINDIF.domaindiscretization=funch{7}{9};
OCMATINDIF.timesettransformation=funch{7}{10};

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
    OCMATCONT.DOMAINDDATA(ii).numode=domaindata(ii).odedim+1;
    OCMATCONT.DOMAINDDATA(ii).eqcoord=1:OCMATCONT.DOMAINDDATA(ii).numeq;
    OCMATCONT.DOMAINDDATA(ii).odecoord=1:domaindata(ii).odedim;
    OCMATCONT.DOMAINDDATA(ii).aecoord=domaindata(ii).odedim+(1:domaindata(ii).aedim);
end

% model information
OCMATINDIF.statecoordinate=statecoord(ocObj);
OCMATINDIF.statecostatecoordinate=1:2*statenum(ocObj);
OCMATINDIF.parametervalue=parametervalue(ocObj);
OCMATINDIF.autonomous=isautonomous(ocObj);

limSet=cell(1,indifforder);
limitsettype=cell(1,indifforder);
limsetarcarg=cell(1,indifforder);
limsetarcint=cell(1,indifforder);
OCMATINDIF.stableflag=cell(1,indifforder);
H=zeros(1,indifforder);
for ii=1:indifforder
    OCMATINDIF.numarc(ii)=arcnum(ocMP(ii));
    Htmp=hamiltonian(ocObj,ocMP(ii));
    H(ii)=Htmp(1);
    if isocasymptotic(ocMP(ii))
        pathtpe{ii}=pathtype(ocMP(ii));
        limSet{ii}=limitset(ocMP(ii));
        if isequilibrium(limSet{ii})
            limitsettype{ii}='e';
        else
            limitsettype{ii}='l';
        end
        OCMATINDIF.stableflag{ii}=~isempty(strfind(pathtpe{ii},'s'));
        limsetarcarg{ii}=arcargument(limSet{ii});
        limsetarcint{ii}=arcinterval(limSet{ii});
    elseif islimitcycle(ocMP(ii))
        pathtpe{ii}='';
        limSet{ii}=ocMP(ii);
        limitsettype{ii}='l';
    else
        return
    end
end
OCMATINDIF.contindifferencevalue=0;
OCMATINDIF.indifferencevalue=H(2)-H(1);
if ~isempty(indifferencevalue)
    OCMATINDIF.indifferencevaluevector=indifferencevalue-(H(2)-H(1));
    OCMATINDIF.contindifferencevalue=1;
end
OCMATINDIF.pathtype=pathtpe;

dimensioncanonicalsystem=canonicalsystemdimension(ocObj);
numberofodes=dimensioncanonicalsystem;

% test if pure state constraints are defined
OCMATINDIF.stateconstraint=stateconstraint(ocObj);
OCMATINDIF.fixendcoordinate=fixendcoordinate;

%
OCMATINDIF.freeendtime=freeendtime;
OCMATINDIF.fixdistance=fixdistance;
OCMATINDIF.fixinitstate=fixinitstate;
for ii=1:indifforder
    if  fixdistance(ii) && isempty(fixdistancecoordinate{ii})
        OCMATINDIF.fixdistancecoordinate{ii}=OCMATINDIF.statecostatecoordinate;
    else
        OCMATINDIF.fixdistancecoordinate{ii}=fixdistancecoordinate{ii};
    end
end

OCMATINDIF.freeparameter=freeparameter;
OCMATINDIF.simple=simple;

OCMATINDIF.targetvalue=targetvalue;
OCMATINDIF.targetcoordinate=targetcoordinate;
OCMATINDIF.targetvectorcoordinate=targetvectorcoordinate;

OCMATINDIF.targetparameter=parameterindex(ocObj,targetparameter);
OCMATINDIF.targetparametervalue=targetparametervalue;

OCMATINDIF.indifferenceorder=indifforder;

sol=generatesolstruct(ocMP,limSet);
parametercoordinate=length(sol.parameters);

ctre=0;
ctrl=0;
octrajectory2limset=zeros(indifforder,2);
for ii=1:indifforder
    if isequilibrium(limSet{ii})
        ctre=ctre+1;
        octrajectory2limset(ii,:)=[ii ctre];
    else
        if ~isempty(OCMATINDIF.freeparameter)
            ctrl=ctrl+1;
            octrajectory2limset(ii,:)=[ii ctrl];
            if ~OCMATINDIF.fixlccoordinate(ii)
                dxdt=canonicalsystem(ocObj,limSet{ii},[],1);

                % for the phase condition
                OCMATINDIF.velocityvector{ctrl}=dxdt(:,1);
                OCMATINDIF.velocitycoordinate{ctrl}=1:length(dxdt(:,1));
                OCMATINDIF.velocityvector{ctrl}=OCMATINDIF.velocityvector{ctrl}/norm(OCMATINDIF.velocityvector{ctrl});
                OCMATINDIF.initialpoint{ctrl}=dependentvar(limSet{ii});
                OCMATINDIF.initialpoint{ctrl}=OCMATINDIF.initialpoint{ctrl}(:,1);
            else
                depvar=dependentvar(limSet{ii});
                OCMATINDIF.velocityvector{ctrl}=[];
                OCMATINDIF.fixlcvalue(ctrl)=depvar(fixlccoordinate(ii),1);

            end
        end
    end
end
OCMATINDIF.equilibriumcounter=ctre;
OCMATINDIF.limitcyclecounter=ctrl;
OCMATINDIF.limitsettype=limitsettype;
OCMATINDIF.limsetarcarg=limsetarcarg;
OCMATINDIF.limsetarcint=limsetarcint;
OCMATINDIF.octrajectory2limset=octrajectory2limset;


for ii=1:OCMATINDIF.indifferenceorder
    if ~isempty(pathtpe{ii})
        if isequilibrium(limSet{octrajectory2limset(ii,2)})
            maptype='c';
        else
            maptype='d';
        end
        J=jacobian(limSet{octrajectory2limset(ii,2)});
        switch OCMATINDIF.limitsettype{ii}
            case 'e'
                if isempty(asymptoticmatrix{ii})
                    Jred=J;
                    [OCMATINDIF.asymptoticmatrix{ii},numstable,numunstable,numcenter,infoStruct]=asymptoticbc(Jred,pathtpe{ii},maptype,ZeroDeviationTolerance(ii));
                    if ~isempty(OCMATINDIF.freeparameter) && ~OCMATINDIF.simple(ii)
                        OCMATINDIF.numstable{ii}=numstable;
                        OCMATINDIF.numunstable{ii}=numunstable;
                        OCMATINDIF.numcenter{ii}=numcenter;
                    end
                else
                    OCMATINDIF.asymptoticmatrix{ii}=asymptoticmatrix{ii};
                end
                limSetdepvar=dependentvar(limSet{ii});
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
                Jred=J;
                [OCMATINDIF.asymptoticmatrix{ii},numstable,numunstable,numcenter,infoStruct]=asymptoticbc(Jred,pathtpe{ii},maptype,ZeroDeviationTolerance(ii));
                if ~isempty(OCMATINDIF.freeparameter) && ~OCMATINDIF.simple(ii)
                    OCMATINDIF.numstable{ii}=numstable;
                    OCMATINDIF.numunstable{ii}=numunstable;
                    OCMATINDIF.numcenter{ii}=numcenter;
                end
                OCMATINDIF.saddlepoint{ii}=dependentvar(limSet{ii});
                OCMATINDIF.saddlepoint{ii}=OCMATINDIF.saddlepoint{ii}(:,end);
        end
    end
end
depvar=dependentvar(ocMP(1));

OCMATINDIF.distance=zeros(1,indifforder);
for ii=1:length(OCMATINDIF.fixdistance)
    if OCMATINDIF.fixdistance(ii)
        tmp=dependentvar(ocMP(ii));
        OCMATINDIF.distance(ii)=norm(OCMATINDIF.saddlepoint{ii}-tmp(OCMATINDIF.statecostatecoordinate,end));
    end
end

OCMATINDIF.cumsumnumarc=cumsum(OCMATINDIF.numarc);
OCMATINDIF.initcoord=[1 OCMATINDIF.cumsumnumarc(1:end-1)+1];

OCMATINDIF.parametervalue=parametervalue(ocObj);
OCMATINDIF.initialtime=sol.x0;

if ~isempty(OCMATINDIF.freeparameter)
    OCMATINDIF.freeparametercoordinate=parametercoordinate+(1:length(OCMATINDIF.freeparameter));
    sol.parameters=[sol.parameters OCMATINDIF.parametervalue(OCMATINDIF.freeparameter)];
    parametercoordinate=length(sol.parameters);
end
OCMATINDIF.freevectorcoordinate=parametercoordinate+(1:size(freevector,2));
for ii=1:size(freevector,2)
    sol.parameters=[sol.parameters 0];
end
if OCMATINDIF.contindifferencevalue
    sol.parameters=[sol.parameters 0];
end
OCMATINDIF.freevector=freevector;
OCMATINDIF.linearization=J;
OCMATINDIF.limitset=limSet;
OCMATINDIF.startvalue=depvar(OCMATINDIF.statecoordinate,1);
pathname=OCMATINDIF.datapath();
[resultfile,globalvarfile]=OCMATINDIF.saveintermediatefiles();
OCMATINDIF.basicresultfilename=fullocmatfile(pathname,resultfile);
OCMATINDIF.basicglobalvarfilename=fullocmatfile(pathname,globalvarfile);
OCMATINDIF.objectivevaluecalc=0;
OCMATINDIF.autonomous=1;

OCMATCONT.HE.numinitialcondition=numel(targetvalue);
OCMATCONT.HE.numendcondition=size(OCMATINDIF.asymptoticmatrix{1},2)*indifforder;
OCMATCONT.numode=numberofodes;
OCMATCONT.maxnumode=max(numberofodes);
OCMATCONT.codimension=1;

function sol=generatesolstruct(ocMP,limSet)
global OCMATINDIF
OCMATINDIF.solutionindex=[];
sol.x=[];
sol.y=[];
sol.arcarg=[];
sol.arcinterval=[];
sol.parameters=[];
ctr=0;
parametercoordinate=0;
for ii=1:OCMATINDIF.indifferenceorder
    if ~isempty(OCMATINDIF.pathtype{ii})
        ctr=ctr+1;
        arcintv=arcinterval(ocMP(ii));
        arcarg=arcargument(ocMP(ii));

        if ctr>1
            sol.x=[sol.x independentvar(ocMP(ii))+sol.x(end)];
        else
            sol.x=independentvar(ocMP(ii));
        end
        sol.y=[sol.y dependentvar(ocMP(ii))];
        sol.arcarg=[sol.arcarg arcarg];
        sol.arcinterval=[sol.arcinterval arcintv];
        sol.x0(ctr)=0;

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
end

OCMATINDIF.limitcycleindex=zeros(1,length(OCMATINDIF.solutionindex));
if ~isempty(OCMATINDIF.freeparameter)
    ctr2=0;
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
end
arcposition=find(diff(sol.x)==0);
sol.arcposition=[[1 arcposition+1];[arcposition numel(sol.x)]];
sol.solver='';
sol.solverinfo.tangent=[];
sol.solverinfo.coeff=[];
sol.solverinfo.tmesh=[];

