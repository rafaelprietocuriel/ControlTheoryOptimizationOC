function sol=initocmat_LAE_LAE(ocObj,ocAsym,contInfo,opt,varargin)
% INITOCMAT_LAE_LAE initialization for limit point solution continuation
% contInfo is the structure returned at the limit point solution

clear global OCMATCONT OCMATLSC
global OCMATCONT OCMATLSC
sol=[];
targetvalue=[];
freevector=[];
targetcoordinate=[];
objectivevaluecalc=[];
freeendtime=[];
pathtpe=pathtype(ocAsym);
if isempty(ocObj)
    ocmatmsg('oc model is empty.')
    return
end
if isempty(ocAsym)
    ocmatmsg('oc trajectory is empty.')
    return
end
if nargin==3
    opt=defaultocoptions;
end
targetvalueidx=find(strcmpi(varargin,'targetvalue'));
targetcoordinateidx=find(strcmpi(varargin,'targetcoordinate'));
freevectoridx=find(strcmpi(varargin,'freevector'));
objectivevaluecalcidx=find(strcmpi(varargin,'objectivevaluecalc'));
freeendtimeidx=find(strcmpi(varargin,'freeendtime'));
if ~isempty(targetvalueidx)
    targetvalue=varargin{targetvalueidx+1};
end
if ~isempty(freeendtimeidx)
    freeendtime=varargin{freeendtimeidx+1};
end
if ~isempty(objectivevaluecalcidx)
    objectivevaluecalc=varargin{objectivevaluecalcidx+1};
end
if ~isempty(targetcoordinateidx)
    targetcoordinate=varargin{targetcoordinateidx+1};
end
if ~isempty(freevectoridx)
    freevector=varargin{freevectoridx+1};
end
if isempty(freeendtime)
    freeendtime=0;
end
if isempty(objectivevaluecalc)
    objectivevaluecalc=false;
end
limSet=limitset(ocAsym);
if isperiodic(limSet)
    limitsettype='d';
    if strcmp(pathtpe,'s')
        freeendtime=-1;
    end
else
    limitsettype='c';
end

OCMATCONT.OPTIONS.SymDerivative=2;
BVPMethod=getocoptions(opt,'GENERAL','BVPMethod');
ZeroDeviationTolerance=getocoptions(opt,'GENERAL','ZeroDeviationTolerance'); % numeric tolerance for a value to be accepted as zero

OCMATCONT.codimension=2;
OCMATLSC.stateconstraint=0;

%OCMATLSC.freeendtime=[];
% initialize global variable (OCMATCONT) for general continuation process 
OCMATCONT.modelname=modelname(ocObj);
%OCMATCONT.modelfunc=modelspecificfunc(ocObj,'4LimitPathContinuation');
OCMATCONT.modelfunc=modelspecificfunc(ocObj,'4SaddlePathContinuation');

% initialize global variable (OCMATLSC) for specific continuation process 
% initialize functions 
% empty function handles mean either that they have to be computed
% numerically (derivatvies, jacobian, hessian) or that standard functions are
% used (e.g. asymptotic boundary condition, initial condition)
funch=OCMATCONT.modelfunc(); % model specific function handles for saddle path continuation

OCMATLSC.canonicalsystem=funch{1};
OCMATLSC.canonicalsystemjacobian=funch{2}{1};
OCMATLSC.canonicalsystemparameterjacobian=funch{2}{2};
OCMATLSC.canonicalsystemhessian=funch{3}{1};
OCMATLSC.canonicalsystemparameterhessian=funch{3}{2};

% function for the boundary conditions
OCMATLSC.bcinitial=funch{5}{1};
OCMATLSC.bcasymptotic=funch{5}{2};
OCMATLSC.bctransversality=funch{5}{3};

% function for Jacobian
OCMATLSC.bcjacobianinitial=funch{6}{1};
OCMATLSC.bcjacobianasymptotic=funch{6}{2};
OCMATLSC.bcjacobiantransversality=funch{6}{3};

% function describing the hybrid structure of the problem
OCMATLSC.hybridinfo=funch{7}{1};
OCMATLSC.domain=funch{7}{2};
OCMATLSC.guard=funch{7}{3};
OCMATLSC.reset=funch{7}{4};
OCMATLSC.switchtime=funch{7}{5};
OCMATLSC.jacobianguard=funch{7}{7};
OCMATLSC.jacobianreset=funch{7}{8};
OCMATLSC.domaindiscretization=funch{7}{9};
OCMATLSC.timesettransformation=funch{7}{10};

% general function
OCMATLSC.plotcontinuation=funch{11};
OCMATLSC.testadmissibility=funch{12};
OCMATLSC.datapath=funch{20};
OCMATLSC.saveintermediatefiles=funch{21};

hybridinfo=OCMATLSC.hybridinfo();
% retrieve information about the arcs from the specifc model file
% domaindata and discretizationdata are ordered in the sequence the
% arcarguments appear in hybridinfo.arcarg
for ii=1:numel(hybridinfo.arcarg)
    domaindata(ii)=OCMATLSC.domain(hybridinfo.arcarg(ii));
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
scoord=statecoord(ocObj);
cscoord=costatecoord(ocObj);
OCMATLSC.statecoord=scoord(:).';
OCMATLSC.statecostatecoord=[scoord(:).' cscoord(:).'];

OCMATLSC.freevector=freevector;
OCMATLSC.targetvalue=targetvalue;
OCMATLSC.targetcoordinate=targetcoordinate;

depvar=dependentvar(ocAsym);

sol=generatesolstruct(ocAsym,BVPMethod);

sol.parameters=sol.arcinterval(2:end-1);
OCMATLSC.switchtimecoord=1:length(sol.parameters);


% mode and path specific variables
limSet=limitset(ocAsym);
J=linearization(limSet);
OCMATLSC.parametervalue=parametervalue(ocObj);
OCMATLSC.initialtime=sol.x0;
OCMATLSC.startvalue=depvar(scoord,1);

OCMATLSC.inftimetransformation=inftimetransformation(ocAsym);
OCMATLSC.truncationtime=sol.arcinterval(end);
OCMATLSC.linearization=J;
OCMATLSC.saddlepoint=dependentvar(limSet);
OCMATLSC.saddlepoint=OCMATLSC.saddlepoint(:,1);
OCMATLSC.asymptoticmatrix=asymptoticbc(J,pathtpe,limitsettype,ZeroDeviationTolerance);
OCMATLSC.pathtype=pathtpe;
pathname=OCMATLSC.datapath();
[resultfile,globalvarfile]=OCMATLSC.saveintermediatefiles();
OCMATLSC.basicresultfilename=fullocmatfile(pathname,resultfile);
OCMATLSC.basicglobalvarfilename=fullocmatfile(pathname,globalvarfile);
OCMATCONT.storedata=1;

if freeendtime
    OCMATLSC.freeendtimecoord=length(sol.parameters)+1;
    sol.parameters=[sol.parameters sol.arcinterval(end)];
    if freeendtime>0
        depvar=dependentvar(ocAsym);
        OCMATLSC.distance=norm(OCMATLSC.saddlepoint-depvar(:,end));
    end
else
    OCMATLSC.truncationtimecoord=[];
end
OCMATLSC.freevectorcoord=length(sol.parameters)+(1:size(freevector,2));
for ii=1:size(freevector,2)
    sol.parameters=[sol.parameters 0];
end

OCMATCONT.HE.numinitialcondition=[];
OCMATCONT.HE.numendcondition=size(OCMATLSC.asymptoticmatrix,2);

OCMATLSC.objectivevaluecalc=objectivevaluecalc;
if objectivevaluecalc && length(sol.y(:,1))==2*statenum(ocObj)
    o=objectivefunction(ocObj,ocAsym,1);
    sol.y(end+1,:)=[0 cumsum((o(1:end-1)+o(2:end))/2.*diff(time(ocObj,ocAsym,1)))];
end
if objectivevaluecalc
    OCMATLSC.objectivevaluecoord=size(sol.y,1);
end

%compute borders
OCMATCONT.OPTIONS.xyvectorized=1;
OCMATCONT.bvpmethod=ocAsym.solver;

if exist('contInfo','var') && ~isempty(contInfo)
    if isfield(contInfo,'data')
        if isfield(contInfo.data,'phi') && isfield(contInfo.data,'psi')
            w=contInfo.data.phi; % right eigenvector to eigenvalue zero
            v=contInfo.data.psi; % left eigenvector to eigenvalue zero
        elseif isfield(contInfo.data,'DFDX')
            % see govaertsetal2005 (p. 242)
            [Q R E]=qr(full(contInfo.data.DFDX));
            R1=R(1:end-1,1:end-1);
            b=R(1:end-1,end);
            w=E*[(R1\-b);1];
            w=w/norm(w);
            v=Q(:,end);
        end
    elseif isfield(contInfo,'LAE_phi')
            w=contInfo.LAE_phi; % right eigenvector to eigenvalue zero
            v=contInfo.LAE_psi; % left eigenvector to eigenvalue zero
    elseif isfield(contInfo,'phi')
            w=contInfo.phi; % right eigenvector to eigenvalue zero
            v=contInfo.psi; % left eigenvector to eigenvalue zero
    end
else
    ocmatmsg('\nTo initialize the bordered matrix set\n\topt=setocoptions(''OCCONTARG'',''WorkSpace'',1);\n')
    w=[];
    v=[];
end

OCMATLSC.LAE_phi=w(:);
OCMATLSC.LAE_psi=v(:);
OCMATLSC.LAE_switch=1;

OCMATLSC.freeendtime=freeendtime;
OCMATLSC.autonomous=isautonomous(ocObj);

function sol=generatesolstruct(ocTrj,solvername,varargin)

sol.x=independentvar(ocTrj);
sol.y=dependentvar(ocTrj);
sol.arcarg=arcargument(ocTrj);
sol.x0=sol.x(1);
sol.arcinterval=arcinterval(ocTrj);
sol.arcposition=arcposition(ocTrj);
sol.parameters=[];
sol.solver=solvername;
