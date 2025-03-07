function sol=initocmat_MM_IS(ocObj,compTrj,continuationtype,varargin)
%
% initocmat_MM_FTE_IS initialization for indifference solution continuation
% of multimodel solution
%
% SOL=initocmat_MM_FTE_IS(MMOBJ,OCTRJM,CONTINUATIONTYPE,CONTINUATIONINDEX,CONTINUATIONTARGET) the
% continuation process is initialized (started) at a solution OCTRJM (cell
% of mmultipath solutions
%
% We assume that the different solutions are of the same structure of parts
% (stages) and the according models are the same for each corresponding
% part
% EXAMPLE: three parts, indifference order two: if the first solution is a
% solution of the models m_1/par_1, m_2/par_2, m_3/par_3, then the second
% solution is of the same sturcture, same models with same parameter values
%
% We will also allow one part/stage solutions of any indifference order:
% this means that for every solution we allow a different model and or
% different parameter values
%
% FREEPARAMETERVECTOR: since the parameter values or/an its location can be
% different for every part in every indifference solution, we have to
% specify the indices of free parameters for each parametervector.
% Example: two parts, indifference order two: alpha,beta appear in all
% parts at the same location, say 1,2 then the freeparameter is given
% by a cellarray of length four: {[1 2]},{[1 2]},{[1 2]},{[1 2]} the target
% value is of the same structure
%
%two parts, indifference order two: alpha,beta appear in different parts at
%a different location´, say part1: [1 2], part2: [2 1], then the
%freeparameter is given by a cellaray of length four: {[1 2]},{[2
%1]},{[1 2]},{[2 1]}
%
%two parts, indifference order two: in the first part we want to continue
%along alpha, beta, [1 2] in the second part along sigma and gamma [3 4]
%then the freeparameter is given by a cellarray of length four:  
% {[1 2]},{[3 4]},{[1 2]},{[3 4]} the target 
% value is of the same structure
% For a more detailed discussion see INITOCMAT_AE_EP
%
clear global OCMATCONT OCMATINDIF
global OCMATCONT OCMATINDIF
sol=[];

fixendstate=[];
fixinitstate=[];
commoninitstate=[];
freeparameter=[];
freestatevector=[];
freestatecoord=[];
objectivevaluecalc=[];
userbc=[];
variationalcalculation=[];
vfreetime=[];
includevariationalobjectivevalue=[];
targetvalue=[];

if isempty(ocObj)
    ocmatmsg('oc model is empty.')
    return
end
if ~(order(ocObj)==order(compTrj))
    ocmatmsg('Model order and solution order have to be equal.')
    return
end
orderofmodel=order(ocObj);
indifferenceorder=orderofmodel;

% contiuationtype: time, state, parameter
% continuationindex: (time) number of switching/endtime, (state)
% coordinate, (parameter) index or name of parameter
% continuationtarget: value
targetvalueidx=find(strcmpi(varargin,'targetvalue'));
fixinitstateidx=find(strcmpi(varargin,'fixinitstate'));
fixendstateidx=find(strcmpi(varargin,'fixendstate'));
commoninitstateidx=find(strcmpi(varargin,'commoninitstate'));
freeparameteridx=find(strcmpi(varargin,'freeparameter'));
freestatevectoridx=find(strcmpi(varargin,'freestatevector'));
freestatecoordidx=find(strcmpi(varargin,'freestatecoord'));
objectivevaluecalcidx=find(strcmpi(varargin,'objectivevaluecalc'));
userbcidx=find(strcmpi(varargin,'userbc'));
variationalcalculationidx=find(strcmpi(varargin,'variationalcalculation'));
vfreetimeidx=find(strcmpi(varargin,'variationalfreetime'));
includevariationalobjectivevalueidx=find(strcmpi(varargin,'includevariationalobjectivevalue'));
if ~isempty(targetvalueidx)
    targetvalue=varargin{targetvalueidx+1};
end
if ~isempty(fixinitstateidx)
    fixinitstate=varargin{fixinitstateidx+1};
end
if ~isempty(commoninitstateidx)
    commoninitstate=varargin{commoninitstateidx+1};
end
if ~isempty(fixendstateidx)
    fixendstate=varargin{fixendstateidx+1};
end
if ~isempty(freeparameteridx)
    freeparameter=varargin{freeparameteridx+1};
end
if ~isempty(freestatevectoridx)
    freestatevector=varargin{freestatevectoridx+1};
end
if ~isempty(freestatecoordidx)
    freestatecoord=varargin{freestatecoordidx+1};
end
if ~isempty(objectivevaluecalcidx)
    objectivevaluecalc=varargin{objectivevaluecalcidx+1};
end
if ~isempty(userbcidx)
    userbc=varargin{userbcidx+1};
end
if ~isempty(variationalcalculationidx)
    variationalcalculation=varargin{variationalcalculationidx+1};
end
if ~isempty(vfreetimeidx)
    vfreetime=varargin{vfreetimeidx+1};
end
if ~isempty(includevariationalobjectivevalueidx)
    includevariationalobjectivevalue=varargin{includevariationalobjectivevalueidx+1};
end

if isempty(objectivevaluecalc)
    objectivevaluecalc=zeros(1,orderofmodel);
end
if isempty(includevariationalobjectivevalue)
    includevariationalobjectivevalue=zeros(1,orderofmodel);
end
if isempty(variationalcalculation)
    variationalcalculation=zeros(1,orderofmodel);
end
if iscell(compTrj)
    compTrj=composite(compTrj);
end
if isempty(userbc)
    userbc=zeros(1,indifferenceorder);
end
funch=cell(1,orderofmodel);
for ii=1:indifferenceorder
    OCMATCONT.modelfunc{ii}=modelspecificfunc(ocObj(ii),'4FiniteHorizonPathContinuation');
    funch{ii}=OCMATCONT.modelfunc{ii}(); % model specific function handles for saddle path continuation
    OCMATCONT.modelname{ii}=modelname(ocObj(ii));

end
OCMATINDIF.orderofmodel=orderofmodel;
OCMATINDIF.indifferenceorder=indifferenceorder;
OCMATINDIF.userbc=userbc;
OCMATINDIF.includevariationalobjectivevalue=includevariationalobjectivevalue;
OCMATINDIF.objectivevaluecalc=objectivevaluecalc;
OCMATINDIF.commoninitstate=commoninitstate;
OCMATINDIF.variationalcalculation=variationalcalculation;

for ii=1:indifferenceorder
    OCMATINDIF.canonicalsystem{ii}=funch{ii}{1};

    OCMATINDIF.canonicalsystemjacobian{ii}=funch{ii}{2}{1};
    OCMATINDIF.canonicalsystemparameterjacobian{ii}=funch{ii}{2}{2};

    % function for the boundary conditions
    OCMATINDIF.bcinitial{ii}=funch{ii}{5}{1};
    OCMATINDIF.bctransversality{ii}=funch{ii}{5}{2};
    OCMATINDIF.bcoptimalhorizon{ii}=funch{ii}{5}{3};
    OCMATINDIF.hitvaluefunc{ii}=funch{ii}{5}{4};
    OCMATINDIF.salvagevalue{ii}=funch{ii}{5}{6};

    OCMATINDIF.guard{ii}=funch{ii}{7}{3};
    OCMATINDIF.reset{ii}=funch{ii}{7}{4};

    if OCMATINDIF.objectivevaluecalc(ii)
        OCMATINDIF.objectivefunction{ii}=funch{ii}{8}{1};
        OCMATINDIF.objectivefunctionjacobian{ii}=funch{ii}{8}{2};
        OCMATINDIF.objectivefunctionparameterjacobian{ii}=funch{ii}{8}{3};
        OCMATINDIF.objectivefunctionderivativetime{ii}=funch{ii}{8}{4};
        if OCMATINDIF.includevariationalobjectivevalue(ii)
            OCMATINDIF.variationalobjectivefunction{ii}=funch{ii}{8}{5};
            OCMATINDIF.variationalobjectivefunctionjacobian{ii}=funch{ii}{8}{6};
            OCMATINDIF.variationalobjectivefunctionparameterjacobian{ii}=funch{ii}{8}{7};
            OCMATINDIF.variationalobjectivefunctionderivativetime{ii}=funch{ii}{8}{8};
            OCMATINDIF.variationalsalvagevalue{ii}=funch{ii}{8}{9};
        end
    end
    if OCMATINDIF.userbc(ii)
        OCMATINDIF.userfunctionbc{ii}=funch{ii}{5}{15};
    end

    if OCMATINDIF.variationalcalculation(ii)
        OCMATINDIF.variationaldynamics{ii}=funch{ii}{4}{4};
        OCMATINDIF.variationaljacobian{ii}=funch{ii}{4}{5};
        OCMATINDIF.variationalparameterjacobian{ii}=funch{ii}{4}{6};
        OCMATINDIF.variationalhamiltonian{ii}=funch{ii}{4}{10};
        OCMATINDIF.variationalguard{ii}=funch{ii}{5}{9};
        OCMATINDIF.variationalreset{ii}=funch{ii}{5}{10};
        OCMATINDIF.variationalbcinitial{ii}=funch{ii}{5}{11};
        OCMATINDIF.variationalbctransversality{ii}=funch{ii}{5}{12};
    end
    %OCMATINDIF.canonicalsystemderivativetime{ii}=funch{ii}{2}{3};
    % general function
    OCMATINDIF.plotcontinuation{ii}=funch{ii}{11};
    OCMATINDIF.testadmissibility{ii}=funch{ii}{12};
    OCMATINDIF.datapath{ii}=funch{ii}{20};
    OCMATINDIF.saveintermediatefiles{ii}=funch{ii}{21};
end
scoord=statecoord(ocObj(1));
cscoord=costatecoord(ocObj(1));
OCMATINDIF.statecostatecoord=[scoord(:).' cscoord(:).'];

sol=generatesolstruct(compTrj);

OCMATINDIF.freestatevector=freestatevector;
arcoffset=0;
switchtimes=cell(1,OCMATINDIF.indifferenceorder);
counter=0;
OCMATINDIF.arcarg=[];
for ii=1:OCMATINDIF.indifferenceorder
    OCMATINDIF.parametervalue{ii}=modelparameter(compTrj(ii));
    arcn=arcnum(compTrj(ii));
    arcintv=arcinterval(compTrj(ii));
    arcarg=arcargument(compTrj(ii));
    arcpos=arcposition(compTrj(ii));
    switchtimes{ii}=arcintv(2:end-1);
    OCMATINDIF.endtime(ii)=arcintv(end);
    OCMATINDIF.edge{ii}=[arcarg(1:end-1);arcarg(2:end)];
    OCMATINDIF.arcarg=[OCMATINDIF.arcarg arcarg];
    OCMATINDIF.numarc(ii)=arcn;
    OCMATINDIF.initialstateindex(ii)=arcpos(end);
    tmp=state(ocObj(ii),compTrj(ii));
    OCMATINDIF.initialstate(:,ii)=tmp(:,1);
    OCMATINDIF.arccoord{ii}=(1:arcn)+arcoffset;
    arcoffset=arcoffset+arcn;
    freeparametercoordinatep1=length(sol.parameters)+1;
    sol.parameters=[sol.parameters switchtimes{ii}];
    OCMATINDIF.switchtimecoord{ii}=freeparametercoordinatep1:length(sol.parameters);
    if variationalcalculation(ii)
        %         if length(vfreetime)~=length(OCMATINDIF.freeswitchingtimeindex)
        %             ocmaterror('Number of variational free time arguments and number of free time arguments are different.')
        %         end
        freeparametercoordinatep1=length(sol.parameters)+1;
        sol.parameters=[sol.parameters vfreetime{ii}];
        OCMATINDIF.vfreetimecoord{ii}=freeparametercoordinatep1:length(sol.parameters);
        OCMATINDIF.freeswitchingtimeindex{ii}=2:2+length(switchtimes{ii})-1;
    end
    counter_start=counter+1;
    counter=counter+OCMATINDIF.numarc(ii);
    OCMATINDIF.solutionindex(counter_start:counter)=ii;
    OCMATINDIF.freemodelparameterindex{ii}=parameterindex(ocObj(ii),freeparameter{ii});
    if ~isempty(freeparameter{ii})
        freeparametercoordinatep1=length(sol.parameters)+1;
        sol.parameters=[sol.parameters OCMATINDIF.parametervalue{ii}(OCMATINDIF.freemodelparameterindex{ii})];
        OCMATINDIF.freemodelparametercoord{ii}=freeparametercoordinatep1:length(sol.parameters);
        OCMATINDIF.initialparametervalue{ii}=OCMATINDIF.parametervalue{ii}(OCMATINDIF.freemodelparameterindex{ii});
    end

end
OCMATINDIF.cumsumnumarc=cumsum(OCMATINDIF.numarc);
OCMATINDIF.initcoord=[1 OCMATINDIF.cumsumnumarc(1:end-1)+1];
OCMATINDIF.endcoord=OCMATINDIF.cumsumnumarc(1:end);
OCMATINDIF.initialtime=sol.x0;
OCMATINDIF.targetvalue=targetvalue;

if ~isempty(freestatevector)
    if isempty(freestatecoord)
        freestatecoord=scoord;
    end
    OCMATINDIF.freestatecoord=freestatecoord;
    freeparametercoordinatep1=length(sol.parameters)+1;
    sol.parameters=[sol.parameters zeros(1,size(freestatevector,2))];
    OCMATINDIF.freestatevectorcoord=freeparametercoordinatep1:length(sol.parameters);
end

OCMATINDIF.fixinitstate=fixinitstate;
OCMATINDIF.fixendstate=fixendstate;

numberofodes=0;
for ii=1:OCMATINDIF.indifferenceorder
    OCMATINDIF.statecoord{ii}=numberofodes+scoord;
    numberofodes=OCMATINDIF.statecoord{ii}(end);
    OCMATINDIF.costatecoord{ii}=numberofodes+scoord;
    numberofodes=OCMATINDIF.costatecoord{ii}(end);
    if OCMATINDIF.variationalcalculation(ii)
        OCMATINDIF.variationaldynamicscoordinate{ii}=numberofodes+(1:cscoord(end));
        numberofodes=numberofodes+cscoord(end);
    end
    if OCMATINDIF.objectivevaluecalc(ii)
        OCMATINDIF.objectivevaluecoord(ii)=numberofodes+1;
        numberofodes=numberofodes+1;
    end
    if OCMATINDIF.includevariationalobjectivevalue(ii)
        OCMATINDIF.variationalobjectivevaluecoord(ii)=numberofodes+1;
        numberofodes=numberofodes+1;
    end
end
switch continuationtype
    case 'initialstate'
        OCMATINDIF.continuationtype=0;
        continuationparameter=[];
    case 'parameter'
        % it is assumed that continuationtarget is of the same structure
        % as continuationindex
        if isempty(targetvalue)
        else
        end
        continuationparameter=0;
end

OCMATINDIF.continuationcoordinate=length(sol.parameters)+1;
sol.parameters=[sol.parameters continuationparameter];

% F canonical system, FV variational system, FO objective function, FE
% exogenous functions, X states/costates, V derivative states/costates, E
% exogenous variables
OCMATINDIF.dFDO=[]; % derivative of the canonical system with respect to the objective variable
OCMATINDIF.dFDV=[]; % derivative of the canonical system with respect to the variational variable
OCMATINDIF.dFDE=[]; % derivative of the canonical system with respect to the exogenous variable
OCMATINDIF.dFDPAR=[]; % derivative of the canonical system with respect to the parameters
OCMATINDIF.dFVDX=[]; % derivative of the variatonal dynamics with respect to the dynamic variable
OCMATINDIF.dFVDO=[]; % derivative of the variatonal dynamics with respect to the objective variable
OCMATINDIF.dFVDE=[]; % derivative of the variatonal dynamics with respect to the exogenous variable
OCMATINDIF.dFVDV=[]; % derivative of the variatonal dynamics with respect to the variatonal variable
OCMATINDIF.dFVDPAR=[]; % derivative of the variatonal dynamics with respect to the parameters
OCMATINDIF.dFODX=[]; % derivative of the objective dynamics with respect to the objective variable
OCMATINDIF.dFODO=[]; % derivative of the objective dynamics with respect to the objective variable
OCMATINDIF.dFODV=[]; % derivative of the objective dynamics with respect to the exogenous variable
OCMATINDIF.dFODE=[]; % derivative of the objective dynamics with respect to the exogenous variable
OCMATINDIF.dFODPAR=[]; % derivative of the objective dynamics with respect to the objective variable
OCMATINDIF.dFEDX=[]; % derivative of the exogenous dynamics with respect to the objective variable
OCMATINDIF.dFEDO=[]; % derivative of the exogenous dynamics with respect to the objective variable
OCMATINDIF.dFEDE=[]; % derivative of the exogenous dynamics with respect to the exogenous variable
OCMATINDIF.dFEDV=[]; % derivative of the exogenous dynamics with respect to the exogenous variable
OCMATINDIF.dFEDPAR=[]; % derivative of the exogenous dynamics with respect to the parameters

dimensioncanonicalsystem=cscoord(end);
OCMATINDIF.dFDPAR=zeros(dimensioncanonicalsystem,length(sol.parameters));
OCMATINDIF.dFDX=zeros(dimensioncanonicalsystem);
OCMATINDIF.dFDXcoord1=1:dimensioncanonicalsystem;
OCMATINDIF.dFDXcoord2=1:dimensioncanonicalsystem;
coord1=dimensioncanonicalsystem;

% coord1 ... coordinates of canonical system / variational system /
% objective dynamics / exogenous dynamics
coord2=0;
if variationalcalculation(1)
    % 1:n ... states, n+1:2n ... costates, 2n+1:3n .... derivative states,
    % 3n+1:4n ... derivative costates
    % 4n+1 ... objective value
    % 4n+2:4n+1+r ... exogenous values
    OCMATINDIF.dFDV=zeros(dimensioncanonicalsystem);
    OCMATINDIF.dFVDX=zeros(dimensioncanonicalsystem);
    OCMATINDIF.dFVDXcoord1=coord1+1:coord1+dimensioncanonicalsystem;
    %OCMATINDIF.dFVDXcoord2=coord2+1:2*dimensioncanonicalsystem;
    OCMATINDIF.dFVDV=zeros(dimensioncanonicalsystem);
    OCMATINDIF.dFVDVcoord1=coord1+1:coord1+dimensioncanonicalsystem;
    OCMATINDIF.dFVDVcoord2=coord2+1:coord2+dimensioncanonicalsystem;
    coord1=coord1+dimensioncanonicalsystem;
    if OCMATINDIF.includevariationalobjectivevalue
        OCMATINDIF.dFVDO=zeros(dimensioncanonicalsystem,2);
    elseif objectivevaluecalc
        OCMATINDIF.dFVDO=zeros(dimensioncanonicalsystem,1);
    end
    OCMATINDIF.dFVDPAR=zeros(dimensioncanonicalsystem,length(sol.parameters));
end
if objectivevaluecalc(1)
    if includevariationalobjectivevalue(1)
        numo=2;
    else
        numo=1;
    end
    OCMATINDIF.dFDO=zeros(dimensioncanonicalsystem,numo);
    OCMATINDIF.dFODO=zeros(numo);
    OCMATINDIF.dFODX=zeros(numo,dimensioncanonicalsystem);
    OCMATINDIF.dFODXcoord1=coord1+1:coord1+numo;
    coord1=coord1+numo;
    OCMATINDIF.dFODXcoord2=1:dimensioncanonicalsystem;
    if variationalcalculation
        OCMATINDIF.dFODV=zeros(numo,dimensioncanonicalsystem);
    end
    OCMATINDIF.dFODPAR=zeros(numo,length(sol.parameters));
end
OCMATINDIF.JX=[[OCMATINDIF.dFDX OCMATINDIF.dFDV OCMATINDIF.dFDO OCMATINDIF.dFDE];[OCMATINDIF.dFVDX OCMATINDIF.dFVDV OCMATINDIF.dFVDO OCMATINDIF.dFVDE];[OCMATINDIF.dFODX OCMATINDIF.dFODV OCMATINDIF.dFODO OCMATINDIF.dFODE];[OCMATINDIF.dFEDX OCMATINDIF.dFEDV OCMATINDIF.dFEDO OCMATINDIF.dFEDE]];
OCMATINDIF.Jpar=zeros(size(sol.y,1),length(sol.parameters));
OCMATINDIF.ODEcoord=1:coord1;

pathname=OCMATINDIF.datapath{1}();
[resultfile,globalvarfile]=OCMATINDIF.saveintermediatefiles{1}();
OCMATINDIF.basicresultfilename=fullocmatfile(pathname,resultfile);
OCMATINDIF.basicglobalvarfilename=fullocmatfile(pathname,globalvarfile);

OCMATCONT.HE.numinitialcondition=[];

OCMATCONT.codimension=1;



function sol=generatesolstruct(ocMultiPath)
global OCMATINDIF
sol.x=[];
sol.y=[];
sol.arcarg=[];
sol.arcinterval=[];
timeshift=0;
for ii=1:OCMATINDIF.indifferenceorder
    ocTrj=ocMultiPath(ii);
    if ii>1
        timeshift=timeshift+arcnum(ocMultiPath(ii-1));
    end
    sol.x=[sol.x independentvar(ocTrj(ii))+timeshift];
    sol.y=[sol.y dependentvar(ocTrj(ii))];
    sol.arcarg=[sol.arcarg arcargument(ocTrj(ii))];
    sol.arcinterval=[sol.arcinterval arcinterval(ocTrj(ii))];
end
sol.parameters=[];
sol.x0=sol.arcinterval(1);
sol.solver='';
sol.parameters=[];
sol.solverinfo.tangent=[];
sol.solverinfo.coeff=[];
sol.solverinfo.tmesh=[];
