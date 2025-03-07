function sol=initocmat_ODE_AE(ocObj,ocAsym,contcoordinate,targetvalue,opt,varargin)
%
% INITOCMAT_AE_AE initialization for asymptotic extremal calculation
%
% SOL=INITOCMAT_AE_AE(OCOBJ,OCASYM,CONTCOORDINATE,TARGETVALUE) the
% continuation process is initialized (started) at a solution OCASYM
% (instance of the ocasymptotic class) that (usually) has been calculated
% during a previous continuation. Therefore, neither the 'TruncationTime' ,
% nor the 'PathType' (see INITOCMAT_AE_EP) can be provided. This
% information is taken from the OCASYM object. 
%
% For a more detailed discussion see INITOCMAT_AE_EP
%
clear global OCMATCONT OCMATAE
global OCMATCONT OCMATAE
sol=[];
freehorizon=[];
freevector=[];
distance=[];
userfile='';
pathtpe=pathtype(ocAsym);
if isempty(ocObj)
    ocmatmsg('oc model is empty.')
    return
end
if isempty(ocAsym)
    ocmatmsg('oc trajectory is empty.')
    return
end
if nargin==4
    opt=defaultocoptions;
end
freehorizonidx=find(strcmpi(varargin,'freehorizon'));
userfileidx=find(strcmpi(varargin,'userfile'));
freevectoridx=find(strcmpi(varargin,'freevector'));
if ~isempty(freehorizonidx)
    freehorizon=varargin{freehorizonidx+1};
end
if ~isempty(userfileidx)
    userfile=varargin{userfileidx+1};
    if ischar(userfile)
        userfile=str2func(userfile);
    end
end
distanceidx=find(strcmpi(varargin,'distance'));
if ~isempty(distanceidx)
    distance=varargin{distanceidx+1};
end
if isempty(freehorizon)
    freehorizon=0;
end
if ~isempty(freevectoridx)
    freevector=varargin{freevectoridx+1};
end

targetvalue=targetvalue(:);

ZeroDeviationTolerance=getocoptions(opt,'GENERAL','ZeroDeviationTolerance'); % numeric tolerance for a value to be accepted as zero

% initialize global variable (OCMATCONT) for general continuation process 
OCMATCONT.modelname=modelname(ocObj);
OCMATCONT.modelfunc=modelspecificfunc(ocObj,'4SaddlePathContinuation');

% initialize global variable (OCMATAE) for specific continuation process 
% initialize functions 
% empty function handles mean either that they have to be computed
% numerically (derivatives, jacobian, hessian) or that standard functions are
% used (e.g. asymptotic boundary condition, initial condition)
funch=OCMATCONT.modelfunc(); % model specific function handles for saddle path continuation

OCMATAE.dynamics=funch{1};
OCMATAE.dynamicsjacobian=funch{2}{1};
OCMATAE.dynamicsparameterjacobian=funch{2}{2};
OCMATAE.dynamicshessian=funch{3}{1};
OCMATAE.dynamicsparameterhessian=funch{3}{2};

% function for the boundary conditions
OCMATAE.bcinitial=funch{5}{1};
OCMATAE.bcasymptotic=funch{5}{2};
try
    OCMATAE.bcdistance=funch{5}{3};
catch
    OCMATAE.bcdistance=[];
end
% function for Jacobian
OCMATAE.bcjacobianinitial=funch{6}{1};
OCMATAE.bcjacobianasymptotic=funch{6}{2};

% function describing the hybrid structure of the problem
OCMATAE.hybridinfo=funch{7}{1};
OCMATAE.domain=funch{7}{2};
OCMATAE.domaindiscretization=funch{7}{3};
OCMATAE.timesettransformation=funch{7}{4};

% general function
OCMATAE.plotcontinuation=funch{11};
OCMATAE.testadmissibility=funch{12};
OCMATAE.datapath=funch{20};
OCMATAE.saveintermediatefiles=funch{21};

hybridinfo=OCMATAE.hybridinfo();
% retrieve information about the arcs from the specifc model file
% domaindata and discretizationdata are ordered in the sequence the
% arcarguments appear in hybridinfo.arcarg
for ii=1:numel(hybridinfo.arcarg)
    domaindata(ii)=OCMATAE.domain(hybridinfo.arcarg(ii));
end

sol=generatesolstruct(ocAsym,solver(ocAsym));

for ii=1:numel(domaindata)
    OCMATCONT.DOMAINDDATA(ii).numode=domaindata(ii).odedim;
    OCMATCONT.DOMAINDDATA(ii).numae=domaindata(ii).aedim;
    OCMATCONT.DOMAINDDATA(ii).daeorder=domaindata(ii).daeorder;
    OCMATCONT.DOMAINDDATA(ii).numeq=numel(domaindata(ii).daeorder);%number of equations
    OCMATCONT.DOMAINDDATA(ii).eqcoord=1:OCMATCONT.DOMAINDDATA(ii).numeq;
    OCMATCONT.DOMAINDDATA(ii).odecoord=1:domaindata(ii).odedim;
    OCMATCONT.DOMAINDDATA(ii).aecoord=domaindata(ii).odedim+(1:domaindata(ii).aedim);
end

% mode and path specific variables
limSet=limitset(ocAsym);
J=linearization(limSet);
depvar=dependentvar(ocAsym);
OCMATAE.parametervalue=parametervalue(ocObj);
OCMATAE.initialtime=sol.x0;
OCMATAE.switchtimecoord=1:numel(sol.parameters)-1;
OCMATAE.inftimetransformation=inftimetransformation(ocAsym);
OCMATAE.truncationtime=sol.arcinterval(end);
OCMATAE.linearization=J;
OCMATAE.saddlepoint=dependentvar(limSet);
OCMATAE.asymptoticmatrix=asymptoticbc(J,pathtpe,'c',ZeroDeviationTolerance);
OCMATAE.targetcoordinate=contcoordinate;
OCMATAE.startvalue=depvar(contcoordinate,1);
OCMATAE.continuationvector=targetvalue-OCMATAE.startvalue;
OCMATAE.pathtype=pathtpe;
OCMATAE.freevector=freevector;
OCMATAE.freevectorindex=1:size(freevector,2);

pathname=OCMATAE.datapath();
[resultfile,globalvarfile]=OCMATAE.saveintermediatefiles();
OCMATAE.basicresultfilename=fullocmatfile(pathname,resultfile);
OCMATAE.basicglobalvarfilename=fullocmatfile(pathname,globalvarfile);

OCMATCONT.HE.numinitialcondition=numel(contcoordinate);
OCMATCONT.HE.numendcondition=size(OCMATAE.asymptoticmatrix,2);

OCMATCONT.codimension=1;

OCMATAE.freehorizon=freehorizon;
if freehorizon
    if ~isempty(distance)
        OCMATAE.distance=distance;
    else
        OCMATAE.distance=norm(sol.y(:,end)-OCMATAE.saddlepoint);
    end
    OCMATAE.freehorizoncoord=length(sol.parameters)+1;
    sol.parameters=[sol.parameters sol.arcinterval(end) zeros(1,size(freevector,2)+1)];
    OCMATCONT.HE.numendcondition=OCMATCONT.HE.numendcondition+1;
else
sol.parameters=[sol.parameters zeros(1,size(freevector,2)+1)];
end

if ~isempty(userfile)
    sol=userfile(sol);
end