function sol=initocmat_AE_IS_P(ocObj,ocMP,parindex,varargin)
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

clear global OCMATCONT OCMATINDIF

global OCMATCONT OCMATINDIF
sol=[];

if isimplicit(ocObj)
    OCMATINDIF.implicit=true;
    sol=initocmat_AE_IS_P4implicit(ocObj,ocMP,parindex,varargin{:});
    return
else
    OCMATINDIF.implicit=false;
end

% input argument ocMP is either a cell of ocasymptotics or a multi path object
ocMP=ocmultipath(ocMP);
indifforder=multiplicity(ocMP);
freevector=[];
fixdistance=[];
freeendtime=[];
simple=[];
opt=[];
targetparametervalue=[];
targetparameterindex=[];
freeparameter=[];
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
    opt=defaultocoptions;
end

stopcriterionidx=find(strcmpi(varargin,'stopcriterion'));
targetparametervalueidx=find(strcmpi(varargin,'targetparametervalue'));
targetparameterindexidx=find(strcmpi(varargin,'targetparameterindex'));
freevectoridx=find(strcmpi(varargin,'freevector'));
freeendtimeidx=find(strcmpi(varargin,'freeendtime'));
fixdistanceidx=find(strcmpi(varargin,'fixdistance')); % can be empty if type is initialstate, if empty for type parameter is set to the state coordinates
simpleidx=find(strcmpi(varargin,'simple'));
optionidx=find(strcmpi(varargin,'option'));
freeparameteridx=find(strcmpi(varargin,'freeparameter'));
userbcidx=find(strcmpi(varargin,'userbc'));
userbcvalueidx=find(strcmpi(varargin,'userbcvalue'));
if ~isempty(stopcriterionidx)
    stopcriterion=varargin{stopcriterionidx+1};
end
if ~isempty(freeendtimeidx)
    freeendtime=varargin{freeendtimeidx+1};
end
if ~isempty(freevectoridx)
    freevector=varargin{freevectoridx+1};
end
if ~isempty(userbcidx)
    userbc=varargin{userbcidx+1};
end
if ~isempty(targetparametervalueidx)
    targetparametervalue=varargin{targetparametervalueidx+1};  
end
if ~isempty(freeparameteridx)
    freeparameter=varargin{freeparameteridx+1};  
end
if ~isempty(userbcvalueidx)
    userbcvalue=varargin{userbcvalueidx+1};  
end
if ~isempty(targetparameterindexidx)
    targetparameterindex=varargin{targetparameterindexidx+1};  
        targetparameterindex=parameterindex(ocObj,targetparameterindex);
end
if ~isempty(fixdistanceidx)
    fixdistance=varargin{fixdistanceidx+1};
end
if isempty(freeendtime)
    freeendtime=zeros(1,indifforder);
elseif numel(freeendtime)==1
    freeendtime=repmat(freeendtime,1,indifforder);
end
if ~isempty(optionidx)
    opt=varargin{optionidx+1};
end
if isempty(opt)
    opt=defaultocoptions;
end

if ~isempty(simpleidx)
    simple=varargin{simpleidx+1};
end
if isempty(freeendtime)
    freeendtime=zeros(1,indifforder); % stop if limitpoint occurs
end
if isempty(userbc)
    userbc=0; % stop if limitpoint occurs
end
if isempty(simple)
    simple=0;
end

if ischar(parindex)
    parindex=parameterindex(ocObj,parindex);
end
if isempty(targetparameterindex) && length(parindex)==1
    targetparameterindex=parindex;
end
if numel(freeendtime)==1
    freeendtime=repmat(freeendtime,1,indifforder);
end

ZeroDeviationTolerance=getocoptions(opt,'GENERAL','ZeroDeviationTolerance'); % numeric tolerance for a value to be accepted as zero
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
% function for the boundary conditions
OCMATINDIF.bcinitial=funch{5}{1};
OCMATINDIF.bcasymptotic=funch{5}{2};
OCMATINDIF.bctransversality=funch{5}{3};
OCMATINDIF.equilibrium=funch{5}{4};
OCMATINDIF.bcindifference=funch{5}{5};
if userbc
    OCMATINDIF.userbcfunc=funch{5}{15};
end

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
if length(simple)==1
    simple=repmat(simple,1,indifforder);
end
limsetcounter=1;
limitcyclecounter=0;
equilibriumcounter=0;
octrajectory2limset=zeros(indifforder,2);
for ii=1:indifforder
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
% model information
OCMATINDIF.statecoordinate=statecoord(ocObj);
OCMATINDIF.statecostatecoordinate=1:2*statenum(ocObj);
OCMATINDIF.parametervalue=parametervalue(ocObj);
OCMATINDIF.simple=simple;
OCMATINDIF.freevector=freevector;
depvar=dependentvar(ocMP(1));
OCMATINDIF.startvalue=depvar(OCMATINDIF.statecoordinate,1);
OCMATINDIF.userbc=userbc;

OCMATINDIF.targetparametervalue=targetparametervalue;
OCMATINDIF.targetparameterindex=parameterindex(ocObj,targetparameterindex);
OCMATINDIF.fixdistance=fixdistance;
OCMATINDIF.pathtype=pathtpe;
OCMATINDIF.freeendtime=freeendtime;
OCMATINDIF.indifforder=indifforder;
OCMATINDIF.limsetcounter=limsetcounter;
OCMATINDIF.equilibriumcounter=equilibriumcounter;
OCMATINDIF.limitcyclecounter=limitcyclecounter;
OCMATINDIF.limitsettype=limitsettype;
OCMATINDIF.limsetarcarg=limsetarcarg;
OCMATINDIF.limsetarcint=limsetarcint;
OCMATINDIF.octrajectory2limset=octrajectory2limset;
OCMATINDIF.indifferenceorder=indifforder;


sol=generatesolstruct(ocMP,limSet);


parametercoordinate=length(sol.parameters);
OCMATINDIF.distance=zeros(1,indifforder);
for ii=1:length(OCMATINDIF.fixdistance)
    OCMATINDIF.saddlepoint{ii}=dependentvar(limSet{octrajectory2limset(ii,2)});
    OCMATINDIF.saddlepoint{ii}=OCMATINDIF.saddlepoint{ii}(:,1);
    if OCMATINDIF.fixdistance(ii)
        OCMATINDIF.distance(ii)=norm(OCMATINDIF.saddlepoint{ii}-ocMP(ii).y(:,end));
    end
end
% mode and path specific variables
ctre=0;
ctrl=0;
for ii=1:limsetcounter
    if isequilibrium(limSet{ii})
        ctre=ctre+1;
        limSetdepvar=dependentvar(limSet{ii});
        sol.parameters=[sol.parameters limSetdepvar(:).'];
        OCMATINDIF.equilibriumcoordinate{ctre}=parametercoordinate+(1:length(limSetdepvar));
        parametercoordinate=length(sol.parameters);
    else
        ctrl=ctrl+1;
        OCMATCONT.monodromy=linearization(limSet{ii});
        dxdt=canonicalsystem(ocObj,limSet{ii},[],1);

        % for the phase condition
        OCMATINDIF.velocityvector{ctrl}=dxdt(:,1);
        OCMATINDIF.velocitycoordinate{ctrl}=1:length(dxdt(:,1));
        OCMATINDIF.velocityvector{ctrl}=OCMATINDIF.velocityvector{ctrl}/norm(OCMATINDIF.velocityvector{ctrl});
        OCMATINDIF.initialpoint{ctrl}=dependentvar(limSet{ii});
        OCMATINDIF.initialpoint{ctrl}=OCMATINDIF.initialpoint{ctrl}(:,1);
    end
end
for ii=1:indifforder
    if isequilibrium(limSet{octrajectory2limset(ii,2)})
        maptype='c';
    else
        maptype='d';
    end
    [asymptoticmatrix,OCMATINDIF.numstable{ii},OCMATINDIF.numunstable{ii},OCMATINDIF.numcenter{ii},infoStruct]=asymptoticbc(linearization(limSet{octrajectory2limset(ii,2)}),pathtpe{ii},maptype,ZeroDeviationTolerance);
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
OCMATINDIF.initcoordinate=[1 OCMATINDIF.cumsumnumarc(1:end-1)+1];
OCMATINDIF.parameterindex=parindex;
OCMATINDIF.initialtime=sol.x0;


%OCMATINDIF.fixcoordinate=fixcoordinate;

OCMATINDIF.freevectorcoordinate=parametercoordinate+(1:size(freevector,2));
sol.parameters=[sol.parameters zeros(1,size(freevector,2))];
parametercoordinate=length(sol.parameters);
OCMATINDIF.parametervaluecoordinate=parametercoordinate+(1:length(parindex));
sol.parameters=[sol.parameters OCMATINDIF.parametervalue(parindex)];
parametercoordinate=length(sol.parameters);
OCMATINDIF.freeparameter=parameterindex(ocObj,freeparameter);
if freeparameter
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

OCMATCONT.codimension=1;
pathname=OCMATINDIF.datapath();
[resultfile,globalvarfile]=OCMATINDIF.saveintermediatefiles();
OCMATINDIF.basicresultfilename=fullocmatfile(pathname,resultfile);
OCMATINDIF.basicglobalvarfilename=fullocmatfile(pathname,globalvarfile);
OCMATINDIF.autonomous=isautonomous(ocObj);


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
for ii=1:OCMATINDIF.indifforder
    ctr=ctr+1;
    arcintv=arcinterval(ocMP(ii));
    arcarg=arcargument(ocMP(ii));

    if ii>1
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



function  Q=computeBase(J,stableflag,NSub)

if ~stableflag
    J=-J;
end
[VU, DU] = realeig(J);
% Select first NSub eigenvectors: unstable eigenspace
% Compute orthonormal basis for the eigenspace
VU = VU(:,1:NSub);
[Q,RU] = qr(VU);


