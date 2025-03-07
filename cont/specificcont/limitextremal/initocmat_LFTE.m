function sol=initocmat_LFTE(ocObj,ocTrj,contInfo,varargin)
% INITOCMAT_LFTE initialization for limit point solution continuation
% contInfo is the structure returned at the limit point solution

clear global OCMATCONT OCMATLSC
global OCMATCONT OCMATLSC
sol=[];
targetvalue=[];
fixendstate=[];
fixinitstate=[];
objectivevaluecalc=[];
targettype='';
freevector=[];
freeparametervector=[];
targetcoordinate=[];
jumpargument=[];
entryindex=[];
exogenousfunction=[];

if isempty(ocObj)
    ocmatmsg('oc model is empty.')
    return
end
if isempty(ocTrj)
    ocmatmsg('oc trajectory is empty.')
    return
end
targetvalueidx=find(strcmpi(varargin,'targetvalue'));
targetcoordinateidx=find(strcmpi(varargin,'targetcoordinate'));
objectivevaluecalcidx=find(strcmpi(varargin,'objectivevaluecalc'));
optimalhorizonidx=find(strcmpi(varargin,'optimalhorizon'));
freevectoridx=find(strcmpi(varargin,'freevector'));
freeparametervectoridx=find(strcmpi(varargin,'freeparametervector'));
targettypeidx=find(strcmpi(varargin,'targettype'));
fixendstateidx=find(strcmpi(varargin,'fixendstate'));
fixinitstateidx=find(strcmpi(varargin,'fixinitstate'));
jumpargumentidx=find(strcmpi(varargin,'jumpargument'));
entryindexidx=find(strcmpi(varargin,'entryindex'));
jumpididx=find(strcmpi(varargin,'jumpid'));
exogenousfunctionidx=find(strcmpi(varargin,'exogenousfunction'));
exogenousinitialstatesidx=find(strcmpi(varargin,'exogenousinitialstates'));
if ~isempty(optimalhorizonidx)
    optimalhorizon=varargin{optimalhorizonidx+1};
else
    optimalhorizon=0;
end
if ~isempty(targetvalueidx)
    targetvalue=varargin{targetvalueidx+1};
end
if ~isempty(targetcoordinateidx)
    targetcoordinate=varargin{targetcoordinateidx+1};
end
if ~isempty(fixendstateidx)
    fixendstate=varargin{fixendstateidx+1};
end
if ~isempty(fixinitstateidx)
    fixinitstate=varargin{fixinitstateidx+1};
end
if ~isempty(jumpargumentidx)
    jumpargument=varargin{jumpargumentidx+1};
end
if ~isempty(entryindexidx)
    entryindex=varargin{entryindexidx+1};
end
if ~isempty(jumpididx)
    jumpid=varargin{jumpididx+1};
end
if ~isempty(freevectoridx)
    freevector=varargin{freevectoridx+1};
end
if ~isempty(freeparametervectoridx)
    freeparametervector=varargin{freeparametervectoridx+1};
end
if ~isempty(exogenousfunctionidx)
    exogenousfunction=varargin{exogenousfunctionidx+1};
end
if ~isempty(exogenousinitialstatesidx)
    exogenousinitialstates=varargin{exogenousinitialstatesidx+1};
end
if ~isempty(objectivevaluecalcidx)
    objectivevaluecalc=varargin{objectivevaluecalcidx+1};
end
if ~isempty(targettypeidx)
    targettype=varargin{targettypeidx+1};
end
if isempty(objectivevaluecalc)
    objectivevaluecalc=false;
end
if isempty(targettype)
    targettype=2;
end
OCMATCONT.OPTIONS.SymDerivative=2;

OCMATCONT.codimension=2;

scoord=statecoord(ocObj);
cscoord=costatecoord(ocObj);
OCMATLSC.statecoord=scoord(:).';
OCMATLSC.statecostatecoord=[scoord(:).' cscoord(:).'];
switch lower(targettype)
    case 'endtime'
        OCMATLSC.targettype=2;
    case 'parameter'
        OCMATLSC.targettype=3;
        parindex=varargin{targettypeidx+2};
        if ischar(parindex)
            parindex=parameterindex(ocObj,parindex);
        end
        OCMATLSC.parameterindex=parindex;
    otherwise
        OCMATLSC.targettype=1;
end
OCMATLSC.freevector=freevector;
OCMATLSC.freeparametervector=freeparametervector;
if ~isempty(freeparametervector)
    OCMATLSC.freeparameterindex=parameterindex(ocObj,freeparametervector);
end

OCMATLSC.exogenousfunction=exogenousfunction;
if exogenousfunction
    OCMATLSC.exogenousinitialstates=exogenousinitialstates;
    OCMATLSC.exogenousnumberofstates=length(exogenousinitialstates);
end

% test if pure state constraints are defined
OCMATLSC.stateconstraint=stateconstraint(ocObj)&&~isempty(jumpargument);

OCMATCONT.modelname=modelname(ocObj);
OCMATCONT.modelfunc=modelspecificfunc(ocObj,'4FiniteHorizonPathContinuation');

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
OCMATLSC.bctransversality=funch{5}{2};
OCMATLSC.bcoptimalhorizon=funch{5}{3};
if OCMATLSC.stateconstraint
    OCMATLSC.bcstateconstraint=funch{5}{7};
    OCMATLSC.bctransversalitysc=funch{5}{8};
end

% function for Jacobian
OCMATLSC.bcjacobianinitial=funch{6}{1};
OCMATLSC.bcjacobiantransversality=funch{6}{2};

% function describing the hybrid structure of the problem
OCMATLSC.hybridinfo=funch{7}{1};
OCMATLSC.domain=funch{7}{2};
OCMATLSC.guard=funch{7}{3};
OCMATLSC.reset=funch{7}{4};
OCMATLSC.switchtime=funch{7}{5};
OCMATLSC.jacobianguard=funch{7}{7};
OCMATLSC.jacobianreset=funch{7}{8};
OCMATLSC.domaindiscretization=funch{7}{9};

if objectivevaluecalc
    OCMATLSC.objectivefunction=funch{8}{1};
    OCMATLSC.objectivefunctionjacobian=funch{8}{2};
    OCMATLSC.objectivefunctionparameterjacobian=funch{8}{3};
    OCMATLSC.objectivefunctionderivativetime=funch{8}{4};
    OCMATLSC.salvagevalue=funch{5}{6};
end
if ~isautonomous(ocObj)
    OCMATLSC.canonicalsystemderivativetime=funch{2}{3};
    OCMATLSC.objectivefunctionderivativetime=funch{8}{4};
end
if OCMATLSC.exogenousfunction
    OCMATLSC.exogenousdynamics=funch{4}{1};
    OCMATLSC.exogenousjacobian=funch{4}{2};
    OCMATLSC.exogenousparameterjacobian=funch{4}{3};
end

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
    if objectivevaluecalc
        OCMATCONT.DOMAINDDATA(ii).numode=domaindata(ii).odedim+1;
        OCMATCONT.DOMAINDDATA(ii).numeq=OCMATCONT.DOMAINDDATA(ii).numeq+1;%number of equations
    end
    OCMATCONT.DOMAINDDATA(ii).eqcoord=1:OCMATCONT.DOMAINDDATA(ii).numeq;
    OCMATCONT.DOMAINDDATA(ii).odecoord=1:domaindata(ii).odedim;
    OCMATCONT.DOMAINDDATA(ii).aecoord=domaindata(ii).odedim+(1:domaindata(ii).aedim);
end
OCMATLSC.parametervalue=parametervalue(ocObj);
OCMATLSC.optimalhorizon=optimalhorizon;

depvar=dependentvar(ocTrj);

sol=generatesolstruct(ocTrj,solver(ocTrj));

arctime=sol.arcinterval;

sol.parameters=sol.arcinterval(2:end-1);
OCMATLSC.switchtimecoord=1:length(sol.parameters);
OCMATLSC.freevectorcoord=length(sol.parameters)+(1:size(freevector,2));
for ii=1:size(freevector,2)
    sol.parameters=[sol.parameters 0];
end

OCMATLSC.transversalityconditioncs=0;
if OCMATLSC.stateconstraint
    OCMATLSC.jumpid=zeros(1,length(sol.arcinterval));
    OCMATLSC.jumpargument=jumpargument;
    if length(jumpid)~=length(jumpargument)
        ocmaterror('Number of jump arguments and identifiers are unequal.')
    end
    OCMATLSC.entryindex=zeros(1,length(sol.arcinterval));
    OCMATLSC.entryindex(entryindex)=1:length(jumpargument);
    OCMATLSC.jumpid(entryindex)=jumpid;
    if ~isempty(jumpargument)
        entrytimecoordinate=length(sol.parameters)+1;
        sol.parameters=[sol.parameters jumpargument];
        OCMATLSC.entrytimecoordinate=entrytimecoordinate:length(sol.parameters);
    end
    if entryindex(end)==length(sol.arcinterval)
        OCMATLSC.transversalityconditioncs=1;
    end
end
if optimalhorizon
    sol.parameters=[sol.parameters  arctime(end)];
    OCMATLSC.endtimecoord=length(sol.parameters);
    OCMATLSC.endtime=[];
end
if ~isempty(freeparametervector)
    if length(OCMATLSC.freeparameterindex)==1 || isempty(OCMATLSC.freeparametertarget)
        freeparametervectorcoordinate=length(sol.parameters)+1;
        sol.parameters=[sol.parameters OCMATLSC.parametervalue(OCMATLSC.freeparameterindex)];
        OCMATLSC.freeparametervectorcoordinate=freeparametervectorcoordinate:length(sol.parameters);
        OCMATLSC.freeinitialparametervalue=[];
    else
        OCMATLSC.freeinitialparametervalue=OCMATLSC.parametervalue(OCMATLSC.freeparameterindex);
        OCMATLSC.freeparameterdirection=OCMATLSC.freeparametertarget-OCMATLSC.freeinitialparametervalue;
        sol.parameters=[sol.parameters 0];
        OCMATLSC.freeparametervectorcoordinate=length(sol.parameters);
    end
end
if OCMATLSC.targettype==2 % end time continuation
    sol.parameters=[sol.parameters  arctime(end)];
    OCMATLSC.endtimecoord=length(sol.parameters);
elseif OCMATLSC.targettype==3
    sol.parameters=[sol.parameters OCMATLSC.parametervalue(parindex)];
    OCMATLSC.parametercoord=length(sol.parameters);
    OCMATLSC.endtime=arctime(end);
else
    OCMATLSC.endtime=arctime(end);
end

OCMATLSC.targetvalue=targetvalue;
OCMATLSC.targetcoordinate=targetcoordinate;
OCMATLSC.fixendstatecoord=fixendstate;
OCMATLSC.fixinitstatecoord=fixinitstate;
if ~isempty(fixendstate)
    OCMATLSC.endstate=depvar(fixendstate,end);
    OCMATLSC.initstate=[];
elseif ~isempty(fixinitstate)
    OCMATLSC.initstate=depvar(fixinitstate,1);
    OCMATLSC.endstate=[];
end

% mode and path specific variables
OCMATLSC.initialtime=sol.x0;
OCMATLSC.startvalue=depvar(scoord,1);
pathname=OCMATLSC.datapath();
[resultfile,globalvarfile]=OCMATLSC.saveintermediatefiles();
OCMATLSC.basicresultfilename=fullocmatfile(pathname,resultfile);
OCMATLSC.basicglobalvarfilename=fullocmatfile(pathname,globalvarfile);
OCMATCONT.storedata=1;

OCMATCONT.HE.numinitialcondition=[];
OCMATCONT.HE.numendcondition=[];

OCMATLSC.objectivevaluecalc=objectivevaluecalc;
dimensioncanonicalsystem=canonicalsystemdimension(ocObj);
numberofodes=dimensioncanonicalsystem;
if objectivevaluecalc && length(sol.y(:,1))==2*statenum(ocObj)
    o=objectivefunction(ocObj,ocTrj,1);
    sol.y(end+1,:)=[0 cumsum((o(1:end-1)+o(2:end))/2.*diff(time(ocObj,ocTrj,1)))];
end
if objectivevaluecalc
    numberofodes=numberofodes+1;
    OCMATLSC.objectivevaluecoord=numberofodes;
else
    sol.y=sol.y(OCMATLSC.statecostatecoord,:);
    OCMATLSC.objectivevaluecoord=[];
end
if OCMATLSC.exogenousfunction
    OCMATLSC.exogenousdynamicscoordinate=numberofodes+1:numberofodes+OCMATLSC.exogenousnumberofstates;
else
    OCMATLSC.exogenousdynamicscoordinate=[];
end


OCMATLSC.dFDO=[]; % derivative of the canonical system with respect to the objective variable
OCMATLSC.dFDE=[]; % derivative of the canonical system with respect to the exogenous variable
OCMATLSC.dFODX=[]; % derivative of the objective dynamics with respect to the objective variable
OCMATLSC.dFODO=[]; % derivative of the objective dynamics with respect to the objective variable
OCMATLSC.dFODE=[]; % derivative of the objective dynamics with respect to the exogenous variable
OCMATLSC.dFODPAR=[]; % derivative of the objective dynamics with respect to the objective variable
OCMATLSC.dFEDX=[]; % derivative of the exogenous dynamics with respect to the objective variable
OCMATLSC.dFEDO=[]; % derivative of the exogenous dynamics with respect to the objective variable
OCMATLSC.dFEDE=[]; % derivative of the exogenous dynamics with respect to the exogenous variable
OCMATLSC.dFEDPAR=[]; % derivative of the exogenous dynamics with respect to the exogenous variable

OCMATLSC.dFDPAR=zeros(dimensioncanonicalsystem,length(sol.parameters));
if objectivevaluecalc
    OCMATLSC.dFDO=zeros(dimensioncanonicalsystem,1);
    OCMATLSC.dFODO=0;
    OCMATLSC.dFODX=zeros(1,dimensioncanonicalsystem);
    if exogenousfunction
        OCMATLSC.dFODE=zeros(1,OCMATLSC.exogenousnumberofstates);
    end
    OCMATLSC.dFODPAR=zeros(1,length(sol.parameters));
end
if OCMATLSC.exogenousfunction
    OCMATLSC.dFDE=zeros(dimensioncanonicalsystem,OCMATLSC.exogenousnumberofstates);
    OCMATLSC.dFEDE=zeros(OCMATLSC.exogenousnumberofstates);
    OCMATLSC.dFEDX=zeros(OCMATLSC.exogenousnumberofstates,dimensioncanonicalsystem);
    if OCMATLSC.objectivevaluecalc
        OCMATLSC.dFEDO=zeros(OCMATLSC.exogenousnumberofstates,1);
    end
    OCMATLSC.dFEDPAR=zeros(OCMATLSC.exogenousnumberofstates,length(sol.parameters));
end
OCMATLSC.Jpar=[OCMATLSC.dFDPAR;OCMATLSC.dFODPAR;OCMATLSC.dFEDPAR];
OCMATLSC.ODEcoord=1:size(OCMATLSC.Jpar,1);
 
OCMATLSC.Jext=[[OCMATLSC.dFDO OCMATLSC.dFDE];[OCMATLSC.dFODO OCMATLSC.dFODE];[OCMATLSC.dFEDO OCMATLSC.dFEDE]];

%compute borders
OCMATCONT.OPTIONS.xyvectorized=1;
OCMATCONT.bvpmethod=ocTrj.solver;

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
