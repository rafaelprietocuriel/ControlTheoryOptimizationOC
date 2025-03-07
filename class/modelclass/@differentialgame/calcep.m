function ocEP=calcep(dgObj,varargin)
%
% CALCEP calculates equilibria of an ocmodel.
%
% CALCEP(OCOBJ) calculates the equilibrium for the canonical system of the
% ocmodel OCOBJ using the symbolic toolbox.
%
% CALCEP(OCOBJ,X0) calculates the equilibrium numerically with starting
% values given at X0. If X0 is empty symbolic calculation is applied. For
% the numerical calculation 'FSOLVE' is used by default.
%
% CALCEP(OCOBJ,X0,ARCARG) equilibria are calculated for the canonical
% system corresponding to arc ARCARG. If ARCARG is empty or missing the
% equilibria for every arc are computed.
%
% CALCEP(OCOBJ,X0,ARCARG,OPT) in the ocoption structure OPT options for the
% numerical solver can be provided (category OPT.EQ). In the category
% GENERAL the equation sover can be changed
% opt.GENERAL.EquationSolver='fsolve' or ...
%
% EP = CALCEP(...) the output argument EP is a cell array of 'dynprimitive'
% instances.

depvar0=[]; % initial guess for equilibria
opt=[]; % option for fsolve
numeric=1;
arcarg=[];
equationfile='';
symbolicequationfile='';
jacobianfile='';

if isempty(dgObj)
    return
end
if nargin>=2
    depvar0=varargin{1};
end
if nargin>=3
    arcarg=varargin{2};
end
if nargin>=4
    opt=varargin{3};
end
if isempty(depvar0)
    % symbolic calculations
    numeric=0;
end
if isempty(opt)
    opt=defaultocoptions;
end

if isempty(arcarg)
    arcarg=arcargument(dgObj);
end
equationfileidx=find(strcmpi(varargin,'equationfile'));
symbolicequationfileidx=find(strcmpi(varargin,'symbolicequationfile'));
jacobianfileidx=find(strcmpi(varargin,'jacobianfile'));
if ~isempty(symbolicequationfileidx)
    symbolicequationfile=varargin{symbolicequationfileidx+1};
end
if ~isempty(equationfileidx)
    equationfile=varargin{equationfileidx+1};
end
if ~isempty(jacobianfileidx)
    jacobianfile=varargin{jacobianfileidx+1};
end
if isempty(symbolicequationfile)
    symbolicequationfile='SymbolicCanonicalSystem';
end

if isempty(equationfile)
    equationfile='CanonicalSystem';
end
if isempty(jacobianfile)
    jacobianfile='CanonicalSystemJacobian';
end

MessageDisplay=strcmp(getocoptions(opt,'INIT','MessageDisplay'),'on');
numarc=numel(arcarg);
if MessageDisplay
    ocmatmsg('\nStart searching for equilibria.\n')
end
counter=0;
if numeric
    numcoord=size(depvar0,1);
    numinitpt=size(depvar0,2);
    numep=numarc*numinitpt;
    epstruct=struct('y',cell(1,numep),'arcarg', cell(1,numep),'solverinfo', cell(1,numep));
else
    [symkernel symbolicinfo]=getsymkernel();
    if isempty(symkernel)
        ocmatmsg('Symbolic toolbox is not activated. Start search using the numeric procedure.\nSearch stopped.\n')
        return
    end
    epstruct=struct('y', cell(1,numarc),'arcarg', cell(1,numarc),'solverinfo', cell(1,numarc));
end
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
    switch numeric
        case 1 % numeric calculation
            numeq=canonicalsystemequationnum(dgObj);
            for ii=1:numinitpt
                if numeq==numcoord
                    epcounter=epcounter+1;
                    [depvar,fval,exitflag,info,jacob]=ocmateqsolve(dgObj,'EquilibriumEquation',0,depvar0(:,ii),act_arcarg,opt);
                    epstruct(epcounter).y=depvar;
                    epstruct(epcounter).arcarg=act_arcarg;
                    epstruct(epcounter).solverinfo=info;
                    epstruct(epcounter).solverinfo.exitflag=exitflag;
                    epstruct(epcounter).solverinfo.fval=fval;
                    epstruct(epcounter).solverinfo.jacobian=jacob;
                else
                    ocmatmsg('Size of initial point ''%d'' is different to number of equations ''%d'' for arc ''%d''.\nCalculation skipped.\n',numcoord,numeq,act_arcarg)
                end
            end

        case 0 % symbolic calculation
            try
                zerofunc=subsparametervalue(dgObj,equilibriumequation(dgObj,[],act_arcarg));
                for ii=1:length(zerofunc)
                    zerofunc(ii)=simplify(zerofunc(ii));
                end
                dependentvariablename=dependentvariable(dgObj,act_arcarg);
                numdependentvar=numel(dependentvariablename);
                sol=ocmatsolve(prepareequationstring(char(zerofunc),symkernel), ...
                    preparevariablestring(cell2vectorstring(dependentvariablename),symkernel),symkernel);
            catch
                sol=[];
            end
            if ~isempty(sol)
                for ii=1:numel(sol)
                    epcounter=epcounter+1;
                    counterdepvar=0;
                    failflag=0;
                    while 1
                        counterdepvar=counterdepvar+1;
                        if counterdepvar>numdependentvar
                            break
                        end
                        try
                            epstruct(epcounter).y(counterdepvar,1)=double(sym(sol(ii).(dependentvariablename{counterdepvar})));
                        catch
                            ocmatmsg('Solution value ''%s'' for ''%s'' is not a numeric value.\nSolution discarded.\n',sol(ii).(dependentvariablename{counterdepvar}),dependentvariablename{counterdepvar})
                            epcounter=epcounter-1;
                            failflag=1;
                            break
                        end
                    end
                    if ~failflag
                        epstruct(epcounter).arcarg=act_arcarg;
                        epstruct(epcounter).solverinfo.Version=symbolicinfo.Version;
                        epstruct(epcounter).solverinfo.solver=symbolicinfo.Name;
                    end
                end
            else
                ocmatmsg('Cannot find solution for arcargument: %d\n',act_arcarg)
            end
    end
end

if numeric % remove empty structure array elements
    epstruct(epcounter+1:numep)=[];
end

% calculate linearization
ocEP=cell(epcounter,1);
for ii=1:epcounter
    ocEP{ii}=dynprimitive(epstruct(ii),dgObj);
    ocEP{ii}=octrajectory(ocEP{ii},linearize(ocEP{ii},dgObj,'dependentvar',1,'jacobianfile',jacobianfile));
end

function eqn=prepareequationstring(eqn,symkernel)

switch symkernel
    case 'maple'
        eqn=['{' removesquarebracket(removematrixstring(eqn)) '}'];
    case 'mupad'
        eqn=['[' removesquarebracket(removematrixstring(eqn)) ']'];
end


function varn=preparevariablestring(varn,symkernel)

switch symkernel
    case 'maple'
        varn=['{' removesquarebracket(varn) '}'];
end