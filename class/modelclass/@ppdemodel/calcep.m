function ocEP=calcep(ppdeObj,depvar0,varargin)
%
% CALCEP calculates an equilibrium solution of a ppdemodel.
%
% CALCEP(PPDEOBJ,X0) calculates the equilibrium numerically with starting
% values given at X0. For the numerical calculation 'FSOLVE' is used by default.
%
% CALCEP(PPDEOBJ,X0,'ARCARG',ARCARG) equilibria are calculated for the canonical
% system corresponding to arc ARCARG. If ARCARG is empty or missing the
% equilibria for every arc are computed.
%
% CALCEP(PPDEOBJ,X0,'OPTION',OPT) in the ocoption structure OPT options for the
% numerical solver can be provided (category OPT.EQ). In the category
% GENERAL the equation sover can be changed
% opt.GENERAL.EquationSolver='fsolve' or ...
%
% EP = CALCEP(...) the output argument EP is a cell array of 'dynprimitive'
% instances.

opt=[]; % option for fsolve
arcarg=[];
equilibriumfile='';
femdata=[];
if isempty(ppdeObj)
    return
end

if nargin==1 || isempty(depvar0)
    fprintf('An initial solution has to be provided. Call\n')
    fprintf('\t\tocEP=calcep(ppdeObj,coeff0,...).\n')
    return
end

arcargidx=find(strcmpi(varargin,'arcarg'));
optidx=find(strcmpi(varargin,'opttion'));
equilibriumfileidx=find(strcmpi(varargin,'equilibriumfile'));
femdataidx=find(strcmpi(varargin,'femdata'));

if ~isempty(arcargidx)
    arcarg=varargin{arcargidx+1};
end
if ~isempty(optidx)
    opt=varargin{optidx+1};
end
if ~isempty(equilibriumfileidx)
    equilibriumfile=varargin{equilibriumfileidx+1};
end
if ~isempty(femdataidx)
    femdata=varargin{femdataidx+1};
end

if isempty(equilibriumfile)
    equilibriumfile='EquilibriumEquation';
end
if isempty(opt)
    opt=defaultocoptions;
end

if isempty(arcarg)
    arcarg=arcargument(ppdeObj);
end
if isempty(femdata)
    n=2*statenum(ppdeObj);
    femdata=generatefemdata(ppdeObj,length(depvar0)/n);
end
MessageDisplay=strcmp(getocoptions(opt,'INIT','MessageDisplay'),'on');
numarc=numel(arcarg);
if MessageDisplay
    ocmatmsg('\nStart searching for equilibria.\n')
end
counter=0;
numcoord=size(depvar0,1);
numinitpt=size(depvar0,2);
numep=numarc*numinitpt;
epstruct=struct('y',cell(1,numep),'arcarg', cell(1,numep),'solverinfo', cell(1,numep));

epcounter=0;
while 1
    counter=counter+1;
    if counter>numarc
        break
    end
    act_arcarg=arcarg(counter);
    if MessageDisplay
        ocmatmsg('\nActual arc identifier for calculation %d\n',act_arcarg)
    end
    for ii=1:numinitpt
        epcounter=epcounter+1;
        [depvar,fval,exitflag,info,jacob]=ocmateqsolve(ppdeObj,equilibriumfile,0,depvar0(:,ii),act_arcarg,opt,femdata);
        epstruct(epcounter).y=depvar;
        epstruct(epcounter).arcarg=act_arcarg;
        epstruct(epcounter).femdata=femdata;
        epstruct(epcounter).solverinfo=info;
        epstruct(epcounter).solverinfo.exitflag=exitflag;
        epstruct(epcounter).solverinfo.fval=fval;
        epstruct(epcounter).solverinfo.jacobian=jacob;
    end
end

epstruct(epcounter+1:numep)=[];


% calculate linearization and generate ppdeprimitive objects
ocEP=cell(epcounter,1);
for ii=1:epcounter
    ocEP{ii}=pdeprimitive(epstruct(ii),ppdeObj);
    ocEP{ii}.linearization=linearize(ocEP{ii},ppdeObj,'dependentvar',1);
end
