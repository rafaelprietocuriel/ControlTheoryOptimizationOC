function sol=initocmat_AE_IS_P4implicit(ocObj,ocMP,parindex,varargin)

global OCMATCONT OCMATINDIF
sol=[];

% input argument ocMP is either a cell of ocasymptotics or a multi path object
% input argument ocMP is either a cell of ocasymptotics or a multi path object
ocMP=ocmultipath(ocMP);
indifforder=multiplicity(ocMP);
freevector=[];
fixdistance=[];
freeendtime=[];
option=[];
targetvalue=[];
targetcoordinate=[];
freeparameter=[];
fixinitstate=[];
targetparametervalue=[];
targetparameter=[];
simple=[];

userbc=[];
userbcvalue=[];
pathtpe=cell(1,indifforder);

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
if nargin==4
    option=defaultocoptions;
end

for ii=1:2:length(varargin)
    if ischar(varargin{ii})
        eval([varargin{ii} '=varargin{ii+1};'])
    end
end

if isempty(option)
    option=defaultocoptions;
end

if isempty(freeendtime)
    freeendtime=zeros(1,indifforder); % stop if limitpoint occurs
end

if isempty(fixdistance)
    fixdistance=zeros(1,indifforder); % stop if limitpoint occurs
end
if isempty(userbc)
    userbc=0; % stop if limitpoint occurs
end
if isempty(simple)
    simple=zeros(1,indifforder);
end
parindex=parameterindex(ocObj,parindex);
if isempty(targetparameter) && length(parindex)==1
    targetparameter=parindex;
end
if length(freeendtime)==1
    freeendtime=repmat(freeendtime,1,indifforder);
end
if length(simple)==1
    simple=repmat(simple,1,indifforder);
end

ZeroDeviationTolerance=zeros(1,indifforder);
for ii=1:indifforder
    if iscell(option)
        ZeroDeviationTolerance(ii)=getocoptions(option{ii},'GENERAL','ZeroDeviationTolerance'); % numeric tolerance for a value to be accepted as zero
    else
        ZeroDeviationTolerance(ii)=getocoptions(option,'GENERAL','ZeroDeviationTolerance'); % numeric tolerance for a value to be accepted as zero
    end
end
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
% function for the boundary conditions
OCMATINDIF.bcinitial=funch{5}{1};
OCMATINDIF.bcasymptotic=funch{5}{2};
OCMATINDIF.bctransversality=funch{5}{3};
OCMATINDIF.equilibrium=funch{5}{4};
OCMATINDIF.bcindifference=funch{5}{5};

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
        if ~all(isident)
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
OCMATINDIF.parametervalue=parametervalue(ocObj);
OCMATINDIF.simple=simple;
OCMATINDIF.freevector=freevector;
depvar=dependentvar(ocMP(1));
OCMATINDIF.startvalue=depvar{1}(OCMATINDIF.statecoordinate,1);
OCMATINDIF.userbc=userbc;
OCMATINDIF.fixinitstate=fixinitstate;
OCMATINDIF.targetvalue=targetvalue;
OCMATINDIF.targetcoordinate=targetcoordinate;

OCMATINDIF.targetparametervalue=targetparametervalue;
OCMATINDIF.targetparameter=parameterindex(ocObj,targetparameter);
OCMATINDIF.fixdistance=fixdistance;
OCMATINDIF.pathtype=pathtpe;
OCMATINDIF.freeendtime=freeendtime;
OCMATINDIF.indifforder=indifforder;
OCMATINDIF.indifferenceorder=indifforder;


sol=generatesolstruct(ocMP,limSet);

parametercoordinate=length(sol.parameters);
ctre=0;
ctrl=0;
for ii=1:limsetcounter
    if isequilibrium(limSet{ii})
        ctre=ctre+1;
        tmp=dependentvar(limSet{ii});
        limSetdepvar=tmp{1};
        sol.parameters=[sol.parameters limSetdepvar(:).'];
        OCMATINDIF.equilibriumcoordinate{ctre}=parametercoordinate+(1:length(limSetdepvar));
        parametercoordinate=length(sol.parameters);
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
for ii=1:indifforder
    if isequilibrium(limSet{octrajectory2limset(ii,2)})
        maptype='c';
    else
        maptype='d';
    end
    tmp=jacobian(limSet{octrajectory2limset(ii,2)});
    J=tmp{1};
    [asymptoticmatrix,OCMATINDIF.numstable{ii},OCMATINDIF.numunstable{ii},OCMATINDIF.numcenter{ii},infoStruct]=asymptoticbc(J,pathtpe{ii},maptype,ZeroDeviationTolerance(ii));
    if OCMATINDIF.simple(ii)
        OCMATINDIF.asymptoticmatrix{ii}=asymptoticmatrix;
    else
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
        %[R11,R12,R21,R22]=RicattiCoeff(OCMATINDIF.Q0{ii},linearization(limSet{octrajectory2limset(ii,2)}),OCMATINDIF.subspacedim{ii});
    end
end
OCMATINDIF.cumsumnumarc=cumsum(OCMATINDIF.numarc);

OCMATINDIF.parameterindex=parindex;
OCMATINDIF.initialtime=sol.x0;

OCMATINDIF.freevectorcoordinate=parametercoordinate+(1:size(freevector,2));
sol.parameters=[sol.parameters zeros(1,size(freevector,2))];
parametercoordinate=length(sol.parameters);
OCMATINDIF.parametervaluecoordinate=parametercoordinate+(1:length(parindex));
sol.parameters=[sol.parameters OCMATINDIF.parametervalue(parindex)];
parametercoordinate=length(sol.parameters);
OCMATINDIF.freeparameter=parameterindex(ocObj,freeparameter);
if ~isempty(freeparameter)
    OCMATINDIF.freeparametercoordinate=parametercoordinate+(1:length(parindex));
    sol.parameters=[sol.parameters OCMATINDIF.parametervalue(OCMATINDIF.freeparameter)];
    parametercoordinate=length(sol.parameters);
end
if userbc
    OCMATINDIF.userbccoordinate=parametercoordinate+1;
    sol.parameters=[sol.parameters 0];
    OCMATINDIF.userbcvalue=userbcvalue;
end
OCMATCONT.HE.numinitialcondition=[];
OCMATCONT.HE.numendcondition=[];

OCMATCONT.numode=numberofodes;
OCMATCONT.maxnumode=max(numberofodes);
OCMATCONT.codimension=1;
pathname=OCMATINDIF.datapath();
[resultfile,globalvarfile]=OCMATINDIF.saveintermediatefiles();
OCMATINDIF.basicresultfilename=fullocmatfile(pathname,resultfile);
OCMATINDIF.basicglobalvarfilename=fullocmatfile(pathname,globalvarfile);
OCMATINDIF.autonomous=isautonomous(ocObj);

function sol=generatesolstruct(ocMultiPath,limSet)
global OCMATINDIF
OCMATINDIF.solutionindex=[];
fn={'x','y','arcarg','arcinterval'};
ctr=0;
parametercoordinate=0;
for ii=1:OCMATINDIF.indifforder
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


% function sol=generatesolstruct(ocMultiPath)
% 
% nummult=multiplicity(ocMultiPath);
% sol=odestruct(ocMultiPath(1));
% fn={'x','y','arcarg','arcinterval'};
% maxodenum=size(sol.y,1);
% for ii=2:nummult
%     soltmp=odestruct(ocMultiPath(ii));
%     maxodenum=max([maxodenum size(soltmp.y,1)]);
%     y=sol.y;
%     sol.y=zeros(maxodenum,size(y,2));
%     sol.y(1:size(y,1),:)=y;
%     y=soltmp.y;
%     soltmp.y=zeros(maxodenum,size(y,2));
%     soltmp.y(1:size(y,1),:)=y;
%     soltmp.x=soltmp.x+sol.x(end);
%     soltmp.arcposition=soltmp.arcposition+length(sol.x);
%     for jj=1:length(fn)
%         sol.(fn{jj})=[sol.(fn{jj}) soltmp.(fn{jj})];
%     end
%     sol.arcpositon=[sol.arcposition soltmp.arcposition];
% end



function  Q=computeBase(J,stableflag,NSub)

if ~stableflag
    J=-J;
end
[VU, DU] = realeig(J);
% Select first NSub eigenvectors: unstable eigenspace
% Compute orthonormal basis for the eigenspace
VU = VU(:,1:NSub);
[Q,RU] = qr(VU);

