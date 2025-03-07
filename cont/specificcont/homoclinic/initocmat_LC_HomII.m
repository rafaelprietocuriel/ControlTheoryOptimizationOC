function sol=initocmat_LC_HomII(ocObj,ocHom,fixedcoord,parindex,opt,varargin)
% INITOCMAT_LC_H_P initialization for the continuation of a limit cycle
% with respect to a parameter, starting from a Hopf bifurcation
%
% SOL=INITOCMAT_LC_H_P(OCOBJ,OCEP,PARINDEX)
%
% SOL=INITOCMAT_LC_H_P(OCOBJ,OCEP,PARINDEX,TARGETVALUE)


clear global OCMATCONT OCMATLC
global OCMATCONT OCMATLC
sol=[];
targetparametervalue=[];
linearizationcalc=[];
if isempty(ocObj)
    ocmatmsg('oc model is empty.')
    return
end
if isempty(ocHom)
    ocmatmsg('oc limitcycle is empty.')
    return
end
% if ~isperiodic(ocHom)
%     ocmatmsg('Second argument is not a limitcycle.')
%     return
% end
if nargin==4
    opt=[];
end
if nargin>=5
    targetparametervalue=varargin{1};
end
if isempty(opt)
    opt=defaultocoptions;
end

if ischar(parindex)
    parindex=parameterindex(ocObj,parindex);
end
% initialize global variable (OCMATCONT) for general continuation process 
OCMATCONT.modelname=modelname(ocObj);
OCMATCONT.modelfunc=modelspecificfunc(ocObj,'4LimitCycleContinuation');

% initialize functions 
% empty function handles mean either that they have to be computed
% numerically (derivatvies, jacobian, hessian) or that standard functions are
% used (e.g. asymptotic boundary condition, initial condition)

funch=OCMATCONT.modelfunc(); % model specific function handles for saddle path continuation

OCMATLC.canonicalsystem=funch{1};
OCMATLC.canonicalsystemjacobian=funch{2}{1};
OCMATLC.canonicalsystemparameterjacobian=funch{2}{2};
OCMATLC.canonicalsystemhessian=funch{3}{1};
OCMATLC.canonicalsystemparameterhessian=funch{3}{2};
% function for the boundary conditions
OCMATLC.bcperiodic=funch{5}{1};

% function for Jacobian
OCMATLC.bcjacobianperiodic=funch{6}{1};

% function describing the hybrid structure of the problem
OCMATLC.hybridinfo=funch{7}{1};
OCMATLC.domain=funch{7}{2};
OCMATLC.guard=funch{7}{3};
OCMATLC.reset=funch{7}{4};
OCMATLC.switchtime=funch{7}{5};
OCMATLC.jacobianguard=funch{7}{7};
OCMATLC.jacobianreset=funch{7}{8};
OCMATLC.domaindiscretization=funch{7}{9};
OCMATLC.timesettransformation=funch{7}{10};
if ~isautonomous(ocObj)
    OCMATLC.canonicalsystemderivativetime=funch{2}{3};
end

% general function
OCMATLC.plotcontinuation=funch{11};
OCMATLC.testadmissibility=funch{12};
OCMATLC.datapath=funch{20};
OCMATLC.saveintermediatefiles=funch{21};

hybridinfo=OCMATLC.hybridinfo();
% retrieve information about the arcs from the specifc model file
% domaindata and discretizationdata are ordered in the sequence the
% arcarguments appear in hybridinfo.arcarg
for ii=1:numel(hybridinfo.arcarg)
    domaindata(ii)=OCMATLC.domain(hybridinfo.arcarg(ii));
end


% test if pure state constraints are defined
OCMATLC.stateconstraint=stateconstraint(ocObj);
sol=generatesolstruct(ocHom);
sccounter=0;
if OCMATLC.stateconstraint
    arcarg4sc=arcargwithactivestateconstraint(ocObj);
    for ii=2:length(sol.arcarg)
        if ismember(sol.arcarg(ii),arcarg4sc)
            sccounter=sccounter+1;
            OCMATLC.jumpcostateindex(sccounter)=ii;
        end
    end
    if sccounter
        OCMATLC.jumpcostatecoord=1:sccounter;
    else
        OCMATLC.jumpcostatecoord=[];
        OCMATLC.jumpcostateindex=[];
    end
end
numswitchtimes=length(sol.arcinterval(2:end-1));

for ii=1:numel(domaindata)
    OCMATCONT.DOMAINDDATA(ii).numode=domaindata(ii).odedim;
    OCMATCONT.DOMAINDDATA(ii).numae=domaindata(ii).aedim;
    OCMATCONT.DOMAINDDATA(ii).daeorder=domaindata(ii).daeorder;
    OCMATCONT.DOMAINDDATA(ii).numeq=numel(domaindata(ii).daeorder);%number of equations
    OCMATCONT.DOMAINDDATA(ii).eqcoord=1:OCMATCONT.DOMAINDDATA(ii).numeq;
    OCMATCONT.DOMAINDDATA(ii).odecoord=1:domaindata(ii).odedim;
    OCMATCONT.DOMAINDDATA(ii).aecoord=domaindata(ii).odedim+(1:domaindata(ii).aedim);
end
sol.parameters=[sol.parameters sol.arcinterval(end) parametervalue(ocObj,parindex)];
% mode and path specific variables
depvar=dependentvar(ocHom);
OCMATLC.parametervalue=parametervalue(ocObj);
OCMATLC.autonomous=isautonomous(ocObj);
OCMATLC.initialtime=sol.x0;
OCMATLC.switchtimecoord=sccounter+[1:numswitchtimes];
OCMATLC.periodcoord=sccounter+1+numswitchtimes;
OCMATLC.parametercoord=sccounter+2+numswitchtimes;
OCMATLC.fixedcoord=fixedcoord;
OCMATLC.fixedvalue=depvar(fixedcoord,1);
OCMATLC.varyparameterindex=parindex;
OCMATLC.targetparametervalue=targetparametervalue;
OCMATLC.linearizationcalc=linearizationcalc;

dxdt=canonicalsystem(ocObj,ocHom,[],1);
OCMATLC.velocityvector=dxdt(:,1);
OCMATLC.velocitycoord=1:length(dxdt(:,1));
OCMATLC.velocityvector=OCMATLC.velocityvector/norm(OCMATLC.velocityvector);
OCMATLC.initialpoint=depvar(:,1);

OCMATCONT.codimension=1;
OCMATCONT.numintegralconstraint=0;

pathname=OCMATLC.datapath();
[resultfile,globalvarfile]=OCMATLC.saveintermediatefiles();
OCMATLC.basicresultfilename=fullocmatfile(pathname,resultfile);
OCMATLC.basicglobalvarfilename=fullocmatfile(pathname,globalvarfile);

function sol=generatesolstruct(ocHom,varargin)

sol.x=independentvar(ocHom);
sol.y=dependentvar(ocHom);
sol.arcarg=arcargument(ocHom);
sol.x0=sol.x(1);
sol.arcinterval=arcinterval(ocHom);
sol.arcposition=arcposition(ocHom);
sol.parameters=parameters(ocHom);
if isempty(sol.parameters)
    sol.parameters=sol.arcinterval(2:end);
end
numcontpar=length(continuationparameter(ocHom));
sol.parameters(end-numcontpar+1:end)=[];
sol.solver='';

sol.solverinfo.tangent=tangent(ocHom);
sol.solverinfo.coeff=ocHom.solverinfo.coeff;
sol.solverinfo.tmesh=ocHom.solverinfo.tmesh;
