function sol=initocmat_ODE_P(odeObj,ocTrj,changeparameterindex,varargin)
% INITOCMAT_AEP_EP initialization for asymptotic extremal calculation
% varying a parameter
%
% SOL=INITOCMAT_AEP_EP(odeObj,ocEP,contidx,targetvalue)


clear global OCMATCONT OCMATODEBVP
global OCMATCONT OCMATODEBVP
sol=[];
targetparametervalue=[];
targetparameterindex=[];
fixinitstate=[];
fixendstate=[];
periodicstate=[];
userbc='';
if isempty(odeObj)
    ocmatmsg('oc model is empty.')
    return
end
if nargin>=4
    targetparametervalue=varargin{1};
end

if ischar(changeparameterindex)
    changeparameterindex=parameterindex(odeObj,changeparameterindex);
end

fixinitstateidx=find(strcmpi(varargin,'fixinitstate'));
fixendstateidx=find(strcmpi(varargin,'fixendstate'));
targetparameterindexidx=find(strcmpi(varargin,'targetparameterindex'));
targetparametervalueidx=find(strcmpi(varargin,'targetparametervalue'));
userbcidx=find(strcmpi(varargin,'userbc'));
periodicstateidx=find(strcmpi(varargin,'periodicstate'));
if ~isempty(fixinitstateidx)
    fixinitstate=varargin{fixinitstateidx+1};
end

if ~isempty(fixendstateidx)
    fixendstate=varargin{fixendstateidx+1};
end
if ~isempty(targetparametervalueidx)
    targetparametervalue=varargin{targetparametervalueidx+1};
end
if ~isempty(targetparameterindexidx)
    targetparameterindex=varargin{targetparameterindexidx+1};
end
if isempty(targetparameterindex)
    targetparameterindex=1:length(targetparametervalue);
end
if ~isempty(periodicstateidx)
    periodicstate=varargin{periodicstateidx+1};
end

if ~isempty(userbcidx)
    userbc=varargin{userbcidx+1};
end

if isempty(userbc)
    userbc=0;
end
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
if userbc
    OCMATODEBVP.bc=funch{5};
end

% function for Jacobian
OCMATODEBVP.bcjacobian=funch{6};
OCMATODEBVP.hybridinfo=funch{7}{1};
OCMATODEBVP.domain=funch{7}{2};
OCMATODEBVP.guard=funch{7}{3};
OCMATODEBVP.reset=funch{7}{4};

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
arctimes=arcinterval(ocTrj);
depvar=dependentvar(ocTrj);

sol=generatesolstruct(ocTrj);

numswitchtimes=length(sol.arcinterval(2:end-1));
sol.parameters=sol.arcinterval(2:end-1);

sol.parameters=[sol.parameters par(changeparameterindex)];

OCMATODEBVP.switchtimecoord=1:numswitchtimes;
OCMATODEBVP.changeparametercoord=numswitchtimes+1:numswitchtimes+length(changeparameterindex);
OCMATODEBVP.autonomous=isautonomous(odeObj);
% mode and path specific variables
OCMATODEBVP.parametervalue=par;
OCMATODEBVP.initialtime=sol.x0;
OCMATODEBVP.switchtimecoord=1:numel(sol.arcinterval)-2;
OCMATODEBVP.targetparameterindex=numswitchtimes+targetparameterindex;
OCMATODEBVP.targetparametervalue=targetparametervalue;
OCMATODEBVP.changeparameterindex=changeparameterindex;
OCMATODEBVP.userbc=userbc;

OCMATODEBVP.fixinitstatecoord=fixinitstate;
OCMATODEBVP.fixendstatecoord=fixendstate;
OCMATODEBVP.periodicstatecoord=periodicstate;
OCMATODEBVP.fixinitialstate=depvar(fixinitstate,1);
OCMATODEBVP.fixendstate=depvar(fixendstate,end);
OCMATODEBVP.endtime=sol.arcinterval(end);

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
