function sol=initocmat_ODE_S(odeObj,ocTrj,contcoordinate,targetvalue,varargin)
% INITOCMAT_ODE_CONT initialization for
%
% SOL=INITOCMAT_ODE_CONT(OCOBJ,SOLINIT,TARGETVALUE)


clear global OCMATCONT OCMATODEBVP
global OCMATCONT OCMATODEBVP

opt=[];
conttype=[]; % 0 ... initial state, 1 ... end state
fixinitstate=[];
fixendstate=[];
periodicstate=[];
freeparameter=[];
usertarget=[];
userbc='';

if isempty(odeObj)
    ocmatmsg('ode model is empty.')
    return
end
if isempty(ocTrj)
    ocmatmsg('ode trajectory is empty.')
    return
end

optionidx=find(strcmpi(varargin,'option'));
conttypeidx=find(strcmpi(varargin,'continuationtype'));
fixinitstateidx=find(strcmpi(varargin,'fixinitstate'));
fixendstateidx=find(strcmpi(varargin,'fixendstate'));
periodicstateidx=find(strcmpi(varargin,'periodicstate'));
freeparameteridx=find(strcmpi(varargin,'freeparameter'));
usertargetidx=find(strcmpi(varargin,'usertarget'));
userbcidx=find(strcmpi(varargin,'userbc'));

if ~isempty(optionidx)
    opt=varargin{optionidx+1};
end
if isempty(opt)
    opt=defaultocoptions;
end

if ~isempty(fixinitstateidx)
    fixinitstate=varargin{fixinitstateidx+1};
end

if ~isempty(fixendstateidx)
    fixendstate=varargin{fixendstateidx+1};
end

if ~isempty(periodicstateidx)
    periodicstate=varargin{periodicstateidx+1};
end

if ~isempty(freeparameteridx)
    freeparameter=varargin{freeparameteridx+1};
end

if ~isempty(userbcidx)
    userbc=varargin{userbcidx+1};
end

if ~isempty(usertargetidx)
    usertarget=varargin{usertargetidx+1};
end
if ischar(freeparameter)
    freeparameter=parameterindex(odeObj,freeparameter);
end

if ~isempty(conttypeidx)
    conttype=varargin{conttypeidx+1};
end
if isempty(conttype)
    conttype=0;
end

if isempty(userbc)
    userbc=0;
end

dim=statenum(odeObj);

%
% initialize global variable (OCMATCONT) for general continuation process 
OCMATCONT.modelname=modelname(odeObj);
OCMATCONT.modelfunc=modelspecificfunc(odeObj,'4Continuation');

% initialize functions 
% empty function handles mean either that they have to be computed
% numerically (derivatvies, jacobian, hessian) or that standard functions are
% used (e.g. asymptotic boundary condition, initial condition)

funch=OCMATCONT.modelfunc(); % model specific function handles for saddle path continuation

OCMATODEBVP.dynamics=funch{1};
OCMATODEBVP.jacobian=funch{2}{1};
OCMATODEBVP.parameterjacobian=funch{2}{2};
OCMATODEBVP.derivativeindependent=funch{2}{2};

% function for the boundary conditions
OCMATODEBVP.bc=funch{5};

% function for Jacobian
OCMATODEBVP.bcjacobian=funch{6};
OCMATODEBVP.hybridinfo=funch{7}{1};
OCMATODEBVP.domain=funch{7}{2};
OCMATODEBVP.guard=funch{7}{3};
OCMATODEBVP.reset=funch{7}{4};
if usertarget
    OCMATODEBVP.usertarget=funch{9};
else
    OCMATODEBVP.usertarget=[];
end
% general function
OCMATODEBVP.plotcontinuation=funch{11};
OCMATODEBVP.testadmissibility=funch{12};
OCMATODEBVP.datapath=funch{20};
OCMATODEBVP.saveintermediatefiles=funch{21};

hybridinfo=OCMATODEBVP.hybridinfo();
% retrieve information about the arcs from the specifc model file
% domaindata and discretizationdata are ordered in the sequence the
% arcarguments appear in hybridinfo.arcarg
for ii=1:numel(hybridinfo.arcarg)
    domaindata(ii)=OCMATODEBVP.domain(hybridinfo.arcarg(ii));
end
par=parametervalue(odeObj);

for ii=1:numel(domaindata)
    OCMATCONT.DOMAINDDATA(ii).numode=domaindata(ii).odedim;
    OCMATCONT.DOMAINDDATA(ii).numae=domaindata(ii).aedim;
    OCMATCONT.DOMAINDDATA(ii).daeorder=domaindata(ii).daeorder;
    OCMATCONT.DOMAINDDATA(ii).numeq=numel(domaindata(ii).daeorder);%number of equations
    OCMATCONT.DOMAINDDATA(ii).eqcoord=1:OCMATCONT.DOMAINDDATA(ii).numeq;
    OCMATCONT.DOMAINDDATA(ii).odecoord=1:domaindata(ii).odedim;
    OCMATCONT.DOMAINDDATA(ii).aecoord=domaindata(ii).odedim+(1:domaindata(ii).aedim);
end
%arctimes=arcinterval(ocTrj);
depvar=dependentvar(ocTrj);

sol=generatesolstruct(ocTrj);

numswitchtimes=length(sol.arcinterval(2:end-1));
sol.parameters=sol.arcinterval(2:end-1);
if ~isempty(freeparameter)
    sol.parameters=[sol.parameters par(freeparameter)];
    OCMATODEBVP.changeparametercoord=numswitchtimes+1:numswitchtimes+length(freeparameter);
end
sol.parameters=[sol.parameters 0];
OCMATODEBVP.switchtimecoord=1:numswitchtimes;
OCMATODEBVP.changeparameterindex=freeparameter;
OCMATODEBVP.autonomous=isautonomous(odeObj);
% mode and path specific variables
OCMATODEBVP.parametervalue=parametervalue(odeObj);
OCMATODEBVP.initialtime=sol.x0;
OCMATODEBVP.arcinterval=sol.arcinterval;

OCMATODEBVP.conttype=conttype;
switch conttype
    case 0
        OCMATODEBVP.startvalue=depvar(contcoordinate,1);
    case 1
        OCMATODEBVP.startvalue=depvar(contcoordinate,end);
end
OCMATODEBVP.contcoordinate=contcoordinate;
OCMATODEBVP.continuationvector=targetvalue-OCMATODEBVP.startvalue;

OCMATODEBVP.statenum=dim;

OCMATODEBVP.fixinitstatecoord=fixinitstate;
OCMATODEBVP.fixendstatecoord=fixendstate;
OCMATODEBVP.periodicstatecoord=periodicstate;
OCMATODEBVP.fixinitialstate=depvar(fixinitstate,1);
OCMATODEBVP.fixendstate=depvar(fixendstate,end);
OCMATODEBVP.endtime=sol.arcinterval(end);
OCMATODEBVP.userbc=userbc;

OCMATCONT.codimension=1;

pathname=OCMATODEBVP.datapath();
[resultfile,globalvarfile]=OCMATODEBVP.saveintermediatefiles();
OCMATODEBVP.basicresultfilename=fullocmatfile(pathname,resultfile);
OCMATODEBVP.basicglobalvarfilename=fullocmatfile(pathname,globalvarfile);

function sol=generatesolstruct(ocTrj)

sol.x=independentvar(ocTrj);
sol.y=dependentvar(ocTrj);
sol.arcarg=arcargument(ocTrj);
sol.x0=sol.x(1);
sol.arcinterval=arcinterval(ocTrj);
sol.arcposition=arcposition(ocTrj);
sol.parameters=[];

if isfield(ocTrj.solverinfo,'tangent')
    sol.solverinfo.tangent=ocTrj.solverinfo.tangent;
else
    sol.solverinfo.tangent=[];
end
