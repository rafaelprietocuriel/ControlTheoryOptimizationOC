function sol=initocmat_AE_FT_IS(mmObj,ocMP,parindex,varargin)
%
clear global OCMATCONT OCMATINDIF
global OCMATCONT OCMATINDIF
sol=[];
freevector=[];
hitstatevalue=[];
hitstatecoordinate=[];
if isempty(mmObj)
    ocmatmsg('oc model is empty.')
    return
end
if isempty(ocMP)
    ocmatmsg('oc trajectory is empty.')
%    return
end
ocMP=ocmultipath(ocMP);

nummod=numberofmodels(mmObj);

optimalhorizonidx=find(strcmpi(varargin,'optimalhorizon'));
fixendstateidx=find(strcmpi(varargin,'fixendstate'));
fixinitstateidx=find(strcmpi(varargin,'fixinitstate'));
freevectoridx=find(strcmpi(varargin,'freevector'));
targetparametervalueidx=find(strcmpi(varargin,'targetparametervalue'));
hitstatevalueidx=find(strcmpi(varargin,'hitstatevalue'));
hitstatecoordinateidx=find(strcmpi(varargin,'hitstatecoordinate'));
if ~isempty(targetparametervalueidx)
    targetparametervalue=varargin{targetparametervalueidx+1};
else
    targetparametervalue=[];
end
if ~isempty(optimalhorizonidx)
    optimalhorizon=varargin{optimalhorizonidx+1};
else
    optimalhorizon=zeros(1,nummod);
end
if ~isempty(fixendstateidx)
    fixendstate=varargin{fixendstateidx+1};
else
    fixendstate=cell(1,nummod);
end
if ~isempty(fixinitstateidx)
    fixinitstate=varargin{fixinitstateidx+1};
else
    fixinitstate=cell(1,nummod);
end
if ~isempty(freevectoridx)
    freevector=varargin{freevectoridx+1};
end
if ~isempty(hitstatevalueidx)
    hitstatevalue=varargin{hitstatevalueidx+1};
end
if ~isempty(hitstatecoordinateidx)
    hitstatecoordinate=varargin{hitstatecoordinateidx+1};
end
if ischar(parindex)
    parindex=parameterindex(mmObj,parindex);
end
if isempty(hitstatevalue)
    hitstatevalue=[];
    hitstatecoordinate=[];
end


% initialize global variable (OCMATCONT) for general continuation process 
OCMATCONT.modelname=modelname(mmObj);
OCMATCONT.modelfunc=modelspecificfunc(mmObj,'4FiniteHorizonPathContinuation');

% initialize global variable (OCMATINDIF) for specific continuation process 
% initialize functions 
% empty function handles mean either that they have to be computed
% numerically (derivatives, jacobian, hessian) or that standard functions are
% used (e.g. asymptotic boundary condition, initial condition)
for ii=1:nummod
    funch{ii}=OCMATCONT.modelfunc{ii}(); % model specific function handles for saddle path continuation
end
for ii=1:nummod
    OCMATINDIF.canonicalsystem{ii}=funch{ii}{1};
    OCMATINDIF.canonicalsystemjacobian{ii}=funch{ii}{2}{1};
    OCMATINDIF.canonicalsystemparameterjacobian{ii}=funch{ii}{2}{2};
    OCMATINDIF.canonicalsystemhessian{ii}=funch{ii}{3}{1};
    OCMATINDIF.canonicalsystemparameterhessian{ii}=funch{ii}{3}{2};

    % function for the boundary conditions
    OCMATINDIF.bcinitial{ii}=funch{ii}{5}{1};
    OCMATINDIF.bctransversality{ii}=funch{ii}{5}{2};
    OCMATINDIF.bcoptimalhorizon{ii}=funch{ii}{5}{3};

    % function for Jacobian
    OCMATINDIF.bcjacobianinitial{ii}=funch{ii}{6}{1};
    OCMATINDIF.bcjacobiantransversality{ii}=funch{ii}{6}{2};

    % function describing the hybrid structure of the problem
    OCMATINDIF.hybridinfo{ii}=funch{ii}{7}{1};
    OCMATINDIF.domain{ii}=funch{ii}{7}{2};
    OCMATINDIF.guard{ii}=funch{ii}{7}{3};
    OCMATINDIF.reset{ii}=funch{ii}{7}{4};
    OCMATINDIF.switchtime{ii}=funch{ii}{7}{5};
    OCMATINDIF.jacobianguard{ii}=funch{ii}{7}{7};
    OCMATINDIF.jacobianreset{ii}=funch{ii}{7}{8};
    OCMATINDIF.domaindiscretization{ii}=funch{ii}{7}{9};

    OCMATINDIF.objectivefunction{ii}=funch{ii}{8}{1};
    OCMATINDIF.objectivefunctionjacobian{ii}=funch{ii}{8}{2};
    OCMATINDIF.objectivefunctionparameterjacobian{ii}=funch{ii}{8}{3};
    OCMATINDIF.salvagevalue{ii}=funch{ii}{5}{6};

    % general function
    OCMATINDIF.plotcontinuation{ii}=funch{ii}{11};
    OCMATINDIF.testadmissibility{ii}=funch{ii}{12};
    OCMATINDIF.datapath{ii}=funch{ii}{20};
    OCMATINDIF.saveintermediatefiles{ii}=funch{ii}{21};

    hybridinfo{ii}=OCMATINDIF.hybridinfo{ii}();

    for jj=1:numel(hybridinfo{ii}.arcarg)
        domaindata{ii}(jj)=OCMATINDIF.domain{ii}(hybridinfo{ii}.arcarg(jj));
    end
    for jj=1:numel(domaindata)
        OCMATCONT.DOMAINDDATA{ii}(jj).numode=domaindata{ii}(jj).odedim;
        OCMATCONT.DOMAINDDATA{ii}(jj).numae=domaindata{ii}(jj).aedim;
        OCMATCONT.DOMAINDDATA{ii}(jj).daeorder=domaindata{ii}(jj).daeorder;
        OCMATCONT.DOMAINDDATA{ii}(jj).numeq=numel(domaindata{ii}(jj).daeorder);%number of equations
        OCMATCONT.DOMAINDDATA{ii}(jj).eqcoord=1:OCMATCONT.DOMAINDDATA{ii}(jj).numeq;
        OCMATCONT.DOMAINDDATA{ii}(jj).odecoord=1:domaindata{ii}(jj).odedim;
        OCMATCONT.DOMAINDDATA{ii}(jj).aecoord=domaindata{ii}(jj).odedim+(1:domaindata{ii}(ii).aedim);
    end
end

OCMATINDIF.targetparametervalue=targetparametervalue;
OCMATINDIF.hitstatevalue=hitstatevalue;
OCMATINDIF.hitstatecoordinate=hitstatecoordinate;

sol=generatesolstruct(ocMP);
% mode and path specific variables
OCMATCONT.HE.edge=cell(1,nummod);
OCMATINDIF.optimalhorizon=optimalhorizon;

arcoffset=0;
parameters=[];
truncationtime=[];
for ii=1:nummod
    arcn=arcnum(ocMP(ii));
    arcintv=arcinterval(ocMP(ii));
    switchtimes{ii}=arcintv(2:end-1);
    truncationtime=[truncationtime arcintv(end)];
    arcarg=arcargument(ocMP(ii));
    OCMATINDIF.edge{ii}=[arcarg(1:end-1);arcarg(2:end)];
    depvar=dependentvar(ocMP(ii));
    OCMATINDIF.fixendstatecoord{ii}=fixendstate{ii};
    OCMATINDIF.fixinitstatecoord{ii}=fixinitstate{ii};
    if ~isempty(fixendstate)
        OCMATINDIF.endstate{ii}=depvar(fixendstate{ii},end);
    end
    if ~isempty(fixinitstate)
        OCMATINDIF.initstate{ii}=depvar(fixinitstate{ii},1);
    end
    OCMATINDIF.numarc(ii)=arcn;
    %OCMATINDIF.initialstateindex(ii)=numel(ocMP(ii).x);
    OCMATINDIF.arccoord{ii}=(1:arcn)+arcoffset;
    arcoffset=arcn;
end
%OCMATINDIF.initialstateindex=cumsum(OCMATINDIF.initialstateindex);
%OCMATINDIF.initialstateindex=[0 OCMATINDIF.initialstateindex(end-1)]+1;
OCMATINDIF.solutionindex=zeros(1,sum(OCMATINDIF.numarc));% relate arcindex to indifference solution

counter=0;

for ii=1:nummod
    if optimalhorizon(ii)
        parameters=[parameters truncationtime(ii)];
        OCMATINDIF.optimalhorizoncoord{ii}=length(parameters);
    else
        OCMATINDIF.truncationtime(ii)=truncationtime(ii);
        OCMATINDIF.optimalhorizoncoord{ii}=[];
    end
    counter_start=counter+1;
    counter=counter+OCMATINDIF.numarc(ii);
    OCMATINDIF.solutionindex(counter_start:counter)=ii;
    if ~isempty(switchtimes{ii})
        OCMATINDIF.switchtimecoord{ii}=length(parameters)+(1:length(switchtimes{ii}));
        parameters=[parameters switchtimes{ii}];
    else
        OCMATINDIF.switchtimecoord{ii}=[];
    end
end
depvar=dependentvar(ocMP(1));

OCMATINDIF.cumsumnumarc=cumsum(OCMATINDIF.numarc);
OCMATINDIF.initcoord=[1 OCMATINDIF.cumsumnumarc(1:end-1)+1];
OCMATINDIF.parameterindex=parindex;
OCMATINDIF.statecoord=statecoord(mmObj);
OCMATINDIF.costatecoord=costatecoord(mmObj);

OCMATINDIF.nummod=nummod;
OCMATINDIF.parametervalue=parametervalue(mmObj);
OCMATINDIF.initialtime=sol.x0;
%lefttimeindex=[1 cumsum(OCMATINDIF.numarc(1:end-1)+1)];
%righttimeindex=cumsum(OCMATINDIF.numarc+1);
%OCMATINDIF.truncationtime=truncationtime;
pathname=OCMATINDIF.datapath{1}();
[resultfile,globalvarfile]=OCMATINDIF.saveintermediatefiles{1}();
OCMATINDIF.basicresultfilename=fullocmatfile(pathname,resultfile);
OCMATINDIF.basicglobalvarfilename=fullocmatfile(pathname,globalvarfile);
OCMATINDIF.objectivevaluecalc=0;
OCMATINDIF.autonomous=1;


%OCMATINDIF.switchtimecoord=length(parameters)+(1:length(sol.parameters));
sol.parameters=[parameters sol.parameters];
OCMATINDIF.parametervaluecoord=length(sol.parameters)+(1:length(parindex{1}));

sol.parameters=[sol.parameters OCMATINDIF.parametervalue{1}(parindex{1})];
OCMATINDIF.freevectorcoord=length(sol.parameters)+(1:size(freevector,2));
for ii=1:size(freevector,2)
    sol.parameters=[sol.parameters 0];
end
OCMATINDIF.freevector=freevector;
OCMATINDIF.startvalue=depvar(OCMATINDIF.statecoord{1},1);

OCMATCONT.HE.numinitialcondition=[];
OCMATCONT.HE.numendcondition=[];

OCMATCONT.codimension=1;


function sol=generatesolstruct(ocMultiPath)
global OCMATCONT
nummult=multiplicity(ocMultiPath);
sol.x=independentvar(ocMultiPath(1));
sol.y=dependentvar(ocMultiPath(1));
sol.y=sol.y(1:OCMATCONT.DOMAINDDATA{1}(1).numode,:);
sol.arcarg=arcargument(ocMultiPath(1));
sol.arcinterval=arcinterval(ocMultiPath(1));
%sol.parameters=parameters(ocMultiPath(1));
%numcontpar=length(continuationparameter(ocMultiPath(1)));
%sol.parameters(end-numcontpar+1:end)=[];
%if isempty(sol.parameters)
%    sol.parameters=sol.arcinterval(2:end-1);
%end
x0=initialtime(ocMultiPath(1));
for ii=2:nummult
    sol.x=[sol.x independentvar(ocMultiPath(ii))-x0+sol.x(end)];
    tmpy=dependentvar(ocMultiPath(ii));
    sol.y=[sol.y tmpy(1:OCMATCONT.DOMAINDDATA{1}(1).numode,:);];
    sol.arcarg=[sol.arcarg arcargument(ocMultiPath(ii))];
    actarcinterval=arcinterval(ocMultiPath(ii));
    sol.arcinterval=[sol.arcinterval actarcinterval];
%    freepar=parameters(ocMultiPath(ii));
%    if isempty(freepar)
%        freepar=actarcinterval(2:end-1);
%        numcontpar=0;
%    else
%        numcontpar=length(continuationparameter(ocMultiPath(ii)));

%    end
%    freepar(end-numcontpar+1:end)=[];
%    sol.parameters=[sol.parameters freepar];
end
sol.parameters=[];
sol.x0=x0;
sol.solver='';

sol.solverinfo.tangent=[];
sol.solverinfo.coeff=[];
sol.solverinfo.tmesh=[];