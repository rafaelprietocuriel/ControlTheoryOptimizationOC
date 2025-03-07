function sol=initocmat_AE_HET4implicit(ocObj,ocMP,parindex,fixcoordinate,varargin)
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
% During the initialization two global variables OCMATCONT and OCMATHET
% are initialized. The first variable contains general information for the
% continuation process. The second variable contains problem specific
% information.
%
% SOL=INITOCMAT_AE_IS(OCOBJ,OCMP,TARGETCOORDINATE,TARGETVALUE,OPT)

global OCMATCONT OCMATHET
sol=[];

if isempty(ocObj)
    ocmatmsg('oc model is empty.')
    return
end
ocMP=ocmultipath(ocMP);
hetorder=multiplicity(ocMP);

if isempty(ocMP)
    ocmatmsg('oc trajectory is empty.')
    return
end

% input argument ocMP is either a cell of ocasymptotics or a multi path object

pathtpe=cell(1,hetorder);
fixdistance=[];
freeendtime=[];
simple=[];
option=[];
hopf=[];
targetparametervalue=[];
targetparameter=[];
fixlccoordinate=zeros(1,hetorder);

for ii=1:2:length(varargin)
    if ischar(varargin{ii})
        eval([varargin{ii} '=varargin{ii+1};'])
    end
end

if ischar(parindex)
    parindex=parameterindex(ocObj,parindex);
end
if isempty(option)
    option=defaultocoptions;
end
if isempty(fixdistance)
    fixdistance=[0 0];
end
if isempty(simple)
    simple=0;
end

ZeroDeviationTolerance=zeros(1,hetorder);
AsymptoticBCmethod=cell(1,hetorder);
for ii=1:hetorder
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

OCMATCONT.modelname=modelname(ocObj);
OCMATCONT.modelfunc=modelspecificfunc(ocObj,'4SaddlePathContinuation');

% initialize global variable (OCMATHET) for specific continuation process
% initialize functions
% empty function handles mean either that they have to be computed
% numerically (derivatvies, jacobian, hessian) or that standard functions
% are
% used (e.g. asymptotic boundary condition, initial condition)
funch=OCMATCONT.modelfunc(); % model specific function handles for saddle path continuation

OCMATHET.canonicalsystem=funch{1};
OCMATHET.canonicalsystemjacobian=funch{2}{1};
OCMATHET.canonicalsystemparameterjacobian=funch{2}{2};
OCMATHET.canonicalsystemhessian=funch{3}{1};
OCMATHET.canonicalsystemparameterhessian=funch{3}{2};

OCMATHET.dimplicitcontroldx=funch{4};

% function for the boundary conditions
OCMATHET.bcinitial=funch{5}{1};
OCMATHET.bcasymptotic=funch{5}{2};
OCMATHET.bctransversality=funch{5}{3};
OCMATHET.equilibrium=funch{5}{4};
OCMATHET.bcindifference=funch{5}{5};
OCMATHET.bclimitcycle=funch{5}{6};

% function for Jacobian
OCMATHET.bcjacobianinitial=funch{6}{1};
OCMATHET.bcjacobianasymptotic=funch{6}{2};
OCMATHET.bcjacobiantransversality=funch{6}{3};

% function describing the hybrid structure of the problem
OCMATHET.hybridinfo=funch{7}{1};
OCMATHET.domain=funch{7}{2};
OCMATHET.guard=funch{7}{3};
OCMATHET.reset=funch{7}{4};
%OCMATHET.switchtime=funch{7}{5};
%OCMATHET.jacobianguard=funch{7}{7};
OCMATHET.jacobianreset=funch{7}{8};
OCMATHET.domaindiscretization=funch{7}{9};
OCMATHET.timesettransformation=funch{7}{10};

% general function
OCMATHET.plotcontinuation=funch{11};
OCMATHET.testadmissibility=funch{12};
OCMATHET.datapath=funch{20};
OCMATHET.saveintermediatefiles=funch{21};

hybridinfo=OCMATHET.hybridinfo();
% retrieve information about the arcs from the specifc model file
% domaindata and discretizationdata are ordered in the sequence the
% arcarguments appear in hybridinfo.arcarg
for ii=1:numel(hybridinfo.arcarg)
    domaindata(ii)=OCMATHET.domain(hybridinfo.arcarg(ii));
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

if length(simple)==1
    simple=repmat(simple,1,hetorder);
end
limsetcounter=1;
limitcyclecounter=0;
equilibriumcounter=0;
octrajectory2limset=zeros(hetorder,2);
limSet=cell(hetorder,1);
limitsettype=cell(hetorder,1);
limsetarcarg=cell(hetorder,1);
limsetarcint=cell(hetorder,1);
numberofodes=[];

for ii=1:hetorder
    numberofodes=[numberofodes odenumber(ocMP(ii))];
    pathtpe{ii}=pathtype(ocMP(ii));
    OCMATHET.stableflag{ii}=~isempty(strfind(pathtpe{ii},'s'));
    if ii>1
        jj=1;
        while ~isidentical(limSet{jj},limitset(ocMP(ii)))
            if jj==limsetcounter
                break
            end
            jj=jj+1;
        end
        if jj==limsetcounter
            limsetcounter=limsetcounter+1;
            limSet{limsetcounter}=limitset(ocMP(ii));
            octrajectory2limset(ii,:)=[ii limsetcounter];
        else
            octrajectory2limset(ii,:)=[ii jj];

        end
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
    OCMATHET.implicitcontrolnum(ii)=length(implicitcontrolcoordinate(ocObj,limsetarcarg{limsetcounter}(end)));
end
for ii=1:limsetcounter
    if isequilibrium(limSet{ii})
        equilibriumcounter=equilibriumcounter+1;
        octrajectory2limset(octrajectory2limset(:,2)==ii,3)=equilibriumcounter;
    else
        limitcyclecounter=limitcyclecounter+1;
        octrajectory2limset(octrajectory2limset(:,2)==ii,3)=limitcyclecounter;
    end

end
limSet(limsetcounter+1:end)=[];
limitsettype(limsetcounter+1:end)=[];
limsetarcarg(limsetcounter+1:end)=[];
limsetarcint(limsetcounter+1:end)=[];

% model information
OCMATHET.statecoordinate=statecoord(ocObj);
OCMATHET.statecostatecoordinate=1:2*statenum(ocObj);
OCMATHET.costatecoordinate=statenum(ocObj)+1:2*statenum(ocObj);
OCMATHET.parametervalue=parametervalue(ocObj);
OCMATHET.simple=simple;
OCMATHET.hopf=hopf;

OCMATHET.targetparametervalue=targetparametervalue;
OCMATHET.targetparameterindex=parameterindex(ocObj,targetparameter);
OCMATHET.fixdistance=fixdistance;
OCMATHET.pathtype=pathtpe;
OCMATHET.freeendtime=freeendtime;
OCMATHET.hetorder=hetorder;
OCMATHET.limsetcounter=limsetcounter;
OCMATHET.equilibriumcounter=equilibriumcounter;
OCMATHET.limitcyclecounter=limitcyclecounter;
OCMATHET.limitsettype=limitsettype;
OCMATHET.limsetarcarg=limsetarcarg;
OCMATHET.limsetarcint=limsetarcint;
OCMATHET.octrajectory2limset=octrajectory2limset;


sol=generatesolstruct(ocMP,limSet);
parametercoordinate=length(sol.parameters);
OCMATHET.distance=zeros(1,hetorder);
for ii=1:length(OCMATHET.fixdistance)
    tmp=dependentvar(limSet{octrajectory2limset(ii,2)});
    OCMATHET.saddlepoint{ii}=tmp{1};
    OCMATHET.saddlepoint{ii}=OCMATHET.saddlepoint{ii}(OCMATHET.statecostatecoordinate);
    if OCMATHET.fixdistance(ii)
        tmp=dependentvar(ocMP(ii));
        OCMATHET.distance(ii)=norm(OCMATHET.saddlepoint{ii}-tmp{end}(OCMATHET.statecostatecoordinate,end));
    end
end
% mode and path specific variables
ctre=0;
ctrl=0;
for ii=1:limsetcounter
    if isequilibrium(limSet{ii})
        ctre=ctre+1;
        tmp=dependentvar(limSet{ii});
        limSetdepvar=tmp{1};
        sol.parameters=[sol.parameters limSetdepvar(:).'];
        OCMATHET.equilibriumcoordinate{ctre}=parametercoordinate+(1:length(limSetdepvar));
        parametercoordinate=length(sol.parameters);
    else
        ctrl=ctrl+1;
        OCMATCONT.monodromy=linearization(limSet{ii});
        dxdt=canonicalsystem(limSet{ii},[],1);
        dxdt=dxdt{1};
        % for the phase condition
        OCMATHET.fixlccoordinate=fixlccoordinate;
        if ~fixlccoordinate(ii)
            OCMATHET.velocityvector{ctrl}=dxdt(:,1);
            OCMATHET.velocitycoordinate{ctrl}=1:length(dxdt(:,1));
            OCMATHET.velocityvector{ctrl}=OCMATHET.velocityvector{ctrl}/norm(OCMATHET.velocityvector{ctrl});
            OCMATHET.initialpoint{ctrl}=dependentvar(limSet{ii});
            OCMATHET.initialpoint{ctrl}=OCMATHET.initialpoint{ctrl}(:,1);
        else
            depvar=dependentvar(limSet{ii});
            depvar=depvar{1};
            OCMATHET.velocityvector{ctrl}=[];
            OCMATHET.fixlccoordinate(ctrl)=fixlccoordinate(ii);
            OCMATHET.fixlcvalue(ctrl)=depvar(fixlccoordinate(ii),1);
        end
    end
end
for ii=1:hetorder
    if isequilibrium(limSet{octrajectory2limset(ii,2)})
        maptype='c';
    else
        maptype='d';
    end
    OCMATHET.maptype{ii}=maptype;
    tmp=jacobian(limSet{octrajectory2limset(ii,2)});
    J=tmp{1};
    [asymptoticmatrix,OCMATHET.numstable{ii},OCMATHET.numunstable{ii},OCMATHET.numcenter{ii},infoStruct]=asymptoticbc(J,pathtpe{ii},maptype,ZeroDeviationTolerance(ii),AsymptoticBCMethod{ii});
    if OCMATHET.simple(ii)
        OCMATHET.asymptoticmatrix{ii}=asymptoticmatrix;
    else
        switch pathtpe{ii}
            case 's'
                OCMATHET.subspacedim{ii}=OCMATHET.numstable{ii};
            case {'u','stu'}
                OCMATHET.subspacedim{ii}=OCMATHET.numunstable{ii};
            case {'sc','cs'}
                OCMATHET.subspacedim{ii}=OCMATHET.numstable{ii}+OCMATHET.numcenter{ii};
        end
        OCMATHET.orthspacedim{ii}=OCMATHET.statecostatecoordinate(end)-OCMATHET.subspacedim{ii};
        Y=zeros(OCMATHET.orthspacedim{ii},OCMATHET.subspacedim{ii});
        OCMATHET.Y{ii}=Y;
        OCMATHET.Q0{ii}=infoStruct.Q;
        OCMATHET.Ycoordinate{ii}=reshape(parametercoordinate+(1:OCMATHET.orthspacedim{ii}*OCMATHET.subspacedim{ii}),OCMATHET.orthspacedim{ii},OCMATHET.subspacedim{ii});
        sol.parameters=[sol.parameters Y(:).'];
        parametercoordinate=length(sol.parameters);
        OCMATHET.Id{ii}=eye(OCMATHET.orthspacedim{ii});
        OCMATHET.numY{ii}=numel(Y);
        %[R11,R12,R21,R22]=RicattiCoeff(OCMATHET.Q0{ii},linearization(limSet{octrajectory2limset(ii,2)}),OCMATHET.subspacedim{ii});
    end
end
OCMATHET.cumsumnumarc=cumsum(OCMATHET.numarc);
OCMATHET.initcoord=[1 OCMATHET.cumsumnumarc(1:end-1)+1];
OCMATHET.parameterindex=parindex;
for ii=1:hetorder/2
    depvar1=dependentvar(ocMP(2*ii-1));
    depvar2=dependentvar(ocMP(2*ii));
    %OCMATHET.initialstatedifference(:,ii)=depvar2{1}(OCMATHET.statecostatecoordinate,1)-depvar1{1}(OCMATHET.statecostatecoordinate,1);
    OCMATHET.initialcostatedifference(:,ii)=depvar2{1}(OCMATHET.costatecoordinate,1)-depvar1{1}(OCMATHET.costatecoordinate,1);
    OCMATHET.fixvalue(:,ii)=depvar1{1}(fixcoordinate,1);
end
OCMATHET.initialtime=sol.x0;


OCMATHET.fixcoordinate=fixcoordinate;
OCMATHET.freestatecoordinate=setdiff(OCMATHET.statecoordinate,fixcoordinate);
pathname=OCMATHET.datapath();
[resultfile,globalvarfile]=OCMATHET.saveintermediatefiles();
OCMATHET.basicresultfilename=fullocmatfile(pathname,resultfile);
OCMATHET.basicglobalvarfilename=fullocmatfile(pathname,globalvarfile);
OCMATHET.autonomous=isautonomous(ocObj);

OCMATHET.parametervaluecoord=length(sol.parameters)+(1:length(parindex));

sol.parameters=[sol.parameters OCMATHET.parametervalue(parindex)];

OCMATHET.findinitconnection=0;
if length(parindex)==hetorder/2
    sol.parameters=[sol.parameters 1];
    OCMATHET.findinitconnection=1;
end

OCMATCONT.HE.numinitialcondition=[];
OCMATCONT.HE.numendcondition=[];

OCMATCONT.numode=numberofodes;
OCMATCONT.maxnumode=max(numberofodes);
OCMATCONT.codimension=1;

function sol=generatesolstruct(ocMultiPath,limSet)
global  OCMATHET  
OCMATHET.solutionindex=[];
fn={'x','y','arcarg','arcinterval'};
ctr=0;
parametercoordinate=0;
for ii=1:OCMATHET.hetorder
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


    OCMATHET.switchtime{ctr}=arcintv(2:end-1);
    OCMATHET.truncationtime(ctr)=arcintv(end);
    OCMATHET.numarc(ctr)=length(arcarg);
    OCMATHET.arcarg{ctr}=arcarg;
    OCMATHET.edge{ctr}=[arcarg(1:end-1);arcarg(2:end)];
    OCMATHET.solutionindex=[OCMATHET.solutionindex repmat(ctr,1,OCMATHET.numarc(ctr))];

    sol.parameters=[sol.parameters OCMATHET.switchtime{ctr}];
    OCMATHET.switchtimecoordinate{ctr}=parametercoordinate+(1:length(OCMATHET.switchtime{ctr}));
    parametercoordinate=length(sol.parameters);
    if OCMATHET.freeendtime(ctr)
        sol.parameters=[sol.parameters OCMATHET.truncationtime(ctr)];
        OCMATHET.truncationtimecoordinate(ctr)=parametercoordinate+1;
        parametercoordinate=length(sol.parameters);
    end
end

ctr2=0;
OCMATHET.limitcycleindex=zeros(1,length(OCMATHET.solutionindex));
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
        
        OCMATHET.switchtime{ctr}=arcintv(2:end-1);
        OCMATHET.numarc(ctr)=length(arcarg);
        OCMATHET.arcarg{ctr}=arcarg;
        OCMATHET.edge{ctr}=[arcarg(1:end-1);arcarg(2:end)];
        OCMATHET.solutionindex=[OCMATHET.solutionindex repmat(ctr,1,OCMATHET.numarc(ctr))];
        OCMATHET.limitcycleindex=[OCMATHET.limitcycleindex repmat(ctr2,1,OCMATHET.numarc(ctr))];

        sol.parameters=[sol.parameters OCMATHET.switchtime{ctr}];
        OCMATHET.switchtimecoordinate{ctr}=parametercoordinate+(1:length(OCMATHET.switchtime{ctr}));
        parametercoordinate=length(sol.parameters);
        
        sol.parameters=[sol.parameters arcintv(end)];
        OCMATHET.periodcoordinate(ctr2)=parametercoordinate+1;
        parametercoordinate=length(sol.parameters);
        
    end
end

arcposition=find(diff(sol.x)==0);
sol.arcposition=[[1 arcposition+1];[arcposition numel(sol.x)]];
sol.solver='';
sol.solverinfo.tangent=[];
sol.solverinfo.coeff=[];
sol.solverinfo.tmesh=[];