function sol=initocmat_ODE_T(odeObj,ocTrj,timeindex,varargin)
% INITOCMAT_ODE_CONT initialization for
%
% SOL=INITOCMAT_ODE_CONT(OCOBJ,SOLINIT,TARGETVALUE)


clear global OCMATCONT OCMATODEBVP
global OCMATCONT OCMATODEBVP

targetvalue=[];
fixinitstate=[];
fixendstate=[];
periodicstate=[];

opt=[];
if isempty(odeObj)
    ocmatmsg('ode model is empty.')
    return
end
if isempty(ocTrj)
    ocmatmsg('ode trajectory is empty.')
    return
end

targetvalueidx=find(strcmpi(varargin,'targetvalue'));
optionidx=find(strcmpi(varargin,'option'));
fixinitstateidx=find(strcmpi(varargin,'fixinitstate'));
fixendstateidx=find(strcmpi(varargin,'fixendstate'));
periodicstateidx=find(strcmpi(varargin,'periodicstate'));
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
if ~isempty(targetvalueidx)
    targetvalue=varargin{targetvalueidx+1};
end

if ~isempty(periodicstateidx)
    periodicstate=varargin{periodicstateidx+1};
end

if isempty(fixinitstate) && isempty(fixendstate) && isempty(periodicstate)
    ocmatmsg('No boundary conditions provided.')
    return
end

dim=statenum(odeObj);

if isnumeric(ocTrj)
    if length(ocTrj)==dim
        n=getocoptions(opt,'GENERAL','TrivialArcMeshNum');
        ocTrjStruct.x=linspace(0,1,n);
        ocTrjStruct.y=ocTrj(:);
        ocTrjStruct.y=ocTrjStruct.y(:,ones(1,n));
        ocTrjStruct.arcarg=0;
        ocTrjStruct.arcinterval=[0 0];
        ocTrj=octrajectory(ocTrjStruct);
    else
        ocmatmsg('initialization vector has wron size.')
        return
    end
end

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

numswitchtimes=length(sol.arcinterval(2:end));
sol.parameters=sol.arcinterval(2:end-1);

sol.parameters=[sol.parameters arctimes(timeindex)];
OCMATODEBVP.switchtimecoord=1:numswitchtimes;
OCMATODEBVP.autonomous=isautonomous(odeObj);
% mode and path specific variables
OCMATODEBVP.parametervalue=parametervalue(odeObj);
OCMATODEBVP.initialtime=sol.x0;
OCMATODEBVP.arcinterval=sol.arcinterval;

OCMATODEBVP.varyparameterindex=timeindex;
OCMATODEBVP.targetvalue=targetvalue;

OCMATODEBVP.statenum=dim;

OCMATODEBVP.fixinitstatecoord=fixinitstate;
OCMATODEBVP.fixendstatecoord=fixendstate;
OCMATODEBVP.periodicstatecoord=periodicstate;
OCMATODEBVP.initialstate=depvar(fixinitstate,1);
OCMATODEBVP.endstate=depvar(fixendstate,end);

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
