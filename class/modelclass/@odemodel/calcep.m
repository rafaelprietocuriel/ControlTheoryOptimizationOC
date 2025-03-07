function ocEP=calcep(odeObj,varargin)
%
% CALCEP calculates equilibria of an ocmodel.
%
% CALCEP(ODEOBJ) calculates the equilibrium for the canonical system of the
% ocmodel ODEOBJ using the symbolic toolbox.
%
% CALCEP(ODEOBJ,X0) calculates the equilibrium numerically with starting
% values given at X0. If X0 is empty symbolic calculation is applied. For
% the numerical calculation 'FSOLVE' is used by default.
%
% CALCEP(ODEOBJ,X0,ARCARG) equilibria are calculated for the canonical
% system corresponding to arc ARCARG. If ARCARG is empty or missing the
% equilibria for every arc are computed.
%
% CALCEP(ODEOBJ,X0,ARCARG,OPT) in the ocoption structure OPT options for the
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

if isempty(odeObj)
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
    arcarg=arcargument(odeObj);
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
            %arcInfoStruct=arcinfo(odeObj,act_arcarg);
            %numeq=arcInfoStruct.odedim+arcInfoStruct.aedim;
            act_arcarg=0;
            for ii=1:numinitpt
                epcounter=epcounter+1;
                [depvar,fval,exitflag,info,jacob]=ocmateqsolve(odeObj,'EquilibriumEquation',depvar0(:,ii),act_arcarg,opt);
                epstruct(epcounter).y=depvar;
                epstruct(epcounter).arcarg=act_arcarg;
                epstruct(epcounter).solverinfo=info;
                epstruct(epcounter).solverinfo.exitflag=exitflag;
                epstruct(epcounter).solverinfo.fval=fval;
                epstruct(epcounter).solverinfo.jacobian=jacob;
            end

        case 0 % symbolic calculation
            try
                zerofunc=subsparametervalue(odeObj,equilibriumequation(odeObj,[],act_arcarg,1));

                dependentvariablename=dependentvariable(odeObj,act_arcarg);
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
                            epstruct(epcounter).y(counterdepvar,1)=double(mystr2sym(sol(ii).(dependentvariablename{counterdepvar})));
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
            end
    end
end

if numeric % remove empty structure array elements
    epstruct(epcounter+1:numep)=[];
end

% calculate linearization
ocEP=cell(epcounter,1);
for ii=1:epcounter
    ocEP{ii}=dynprimitive(epstruct(ii),odeObj);
    ocEP{ii}.octrajectory.linearization=equilibriumjacobian(ocEP{ii},odeObj,'dependentvar');
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