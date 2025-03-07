function sol=initocmat_LAE_LAE_LC(ocObj,ocAsym,contcoordinate,fixedcoordinate,opt,contInfo,varargin)
% INITOCMAT_LAE_LAE initialization for limit point solution continuation
% contInfo is the structurre returned at the limit point solution

clear global OCMATCONT OCMATLSC
global OCMATCONT OCMATLSC
sol=[];
targetvalue=[];
pathtpe=pathtype(ocAsym);
if isempty(ocObj)
    ocmatmsg('oc model is empty.')
    return
end
if isempty(ocAsym)
    ocmatmsg('oc trajectory is empty.')
    return
end
if nargin==5
    opt=defaultocoptions;
end
targetvalueidx=find(strcmpi(varargin,'targetvalue'));
if ~isempty(targetvalueidx)
    targetvalue=varargin{targetvalueidx+1};
end
OCMATCONT.OPTIONS.SymDerivative=2;
BVPMethod=getocoptions(opt,'GENERAL','BVPMethod');
ZeroDeviationTolerance=getocoptions(opt,'GENERAL','ZeroDeviationTolerance'); % numeric tolerance for a value to be accepted as zero

if 0%~strcmp(solver(ocAsym),BVPMethod)
    % remove solverinfo since actual BVP solver and solver used for the
    % calculation of ocAsym are not identical
    ocAsym.solver='';
    ocAsym.solverinfo=[];
end
OCMATCONT.codimension=2;
OCMATLSC.stateconstraint=0;
OCMATLSC.objectivevaluecalc=0;

%OCMATLSC.movinghorizon=[];
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

OCMATLSC.freecoordinate=setdiff(1:statenum(ocObj),[contcoordinate fixedcoordinate]);
depvar=dependentvar(ocAsym);
sol=generatesolstruct(ocAsym,BVPMethod);
sol.parameters=[sol.parameters(1:end-1) depvar(OCMATLSC.freecoordinate,1) depvar(contcoordinate,1).'];


% mode and path specific variables
limSet=limitset(ocAsym);
J=linearization(limSet);
OCMATLSC.parametervalue=parametervalue(ocObj);
OCMATLSC.initialtime=sol.x0;
OCMATLSC.switchtimecoord=1:numel(sol.parameters)-2;
OCMATLSC.freecoordinateindex=length(OCMATLSC.switchtimecoord)+1;
OCMATLSC.inftimetransformation=inftimetransformation(ocAsym);
OCMATLSC.truncationtimecoord=length(OCMATLSC.switchtimecoord);
OCMATLSC.linearization=J;
OCMATLSC.saddlepoint=dependentvar(limSet);
OCMATLSC.saddlepoint=OCMATLSC.saddlepoint(:,end);
OCMATLSC.asymptoticmatrix=asymptoticbc(J,pathtpe,'d',ZeroDeviationTolerance);
OCMATLSC.contcoordinate=contcoordinate;
OCMATLSC.targetvalue=targetvalue;
OCMATLSC.fixedcoordinate=fixedcoordinate;
OCMATLSC.startvalue=depvar(sort([contcoordinate OCMATLSC.freecoordinate]),1);
OCMATLSC.fixedstate=depvar(fixedcoordinate,1);
OCMATLSC.pathtype=pathtpe;
pathname=OCMATLSC.datapath();
[resultfile,globalvarfile]=OCMATLSC.saveintermediatefiles();
OCMATLSC.basicresultfilename=fullocmatfile(pathname,resultfile);
OCMATLSC.basicglobalvarfilename=fullocmatfile(pathname,globalvarfile);
OCMATCONT.storedata=1;

OCMATCONT.HE.numinitialcondition=numel(contcoordinate);
OCMATCONT.HE.numendcondition=size(OCMATLSC.asymptoticmatrix,2);

%compute borders
OCMATCONT.OPTIONS.xyvectorized=1;
OCMATCONT.bvpmethod=ocAsym.solver;

if ~isempty(contInfo)
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
    else
    end
else
    ocmaterror('Not implemented yte.')
end

OCMATLSC.LAE_phi=w(:);
OCMATLSC.LAE_psi=v(:);
OCMATLSC.LAE_switch=1;


function sol=generatesolstruct(ocTrj,solvername,varargin)

sol.x=independentvar(ocTrj);
sol.y=dependentvar(ocTrj);
sol.arcarg=arcargument(ocTrj);
sol.x0=sol.x(1);
sol.arcinterval=arcinterval(ocTrj);
sol.arcposition=arcposition(ocTrj);
sol.parameters=parameters(ocTrj);
if isfield(ocTrj.solverinfo,'conttype')
    switch ocTrj.solverinfo.conttype
        case 'extremal2ep'
            if ~isempty(sol.parameters)
                numcontpar=length(continuationparameter(ocTrj));
                sol.parameters(end-numcontpar+1:end)=[];
            end
        otherwise
            sol.parameters=[];
    end
end
sol.solver=solvername;

if isfield(ocTrj.solverinfo,'tangent')
    sol.solverinfo.tangent=ocTrj.solverinfo.tangent;
else
    sol.solverinfo.tangent=[];
end
if isfield(ocTrj.solverinfo,'coeff')
    sol.solverinfo.coeff=ocTrj.solverinfo.coeff;
else
    sol.solverinfo.coeff=[];
end
if isfield(ocTrj.solverinfo,'yp')
    sol.solverinfo.yp=ocTrj.solverinfo.yp;
end
if isfield(ocTrj.solverinfo,'ypmid')
    sol.solverinfo.ypmid=ocTrj.solverinfo.ypmid;
end
if isfield(ocTrj.solverinfo,'ymid')
    sol.solverinfo.yp=ocTrj.solverinfo.yp;
end
if isfield(ocTrj.solverinfo,'ymid')
    sol.solverinfo.ymid=ocTrj.solverinfo.ymid;
end