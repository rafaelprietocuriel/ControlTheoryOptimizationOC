function ocFP=calcfp(docObj,varargin)
%
% CALCFP calculates fix point of a discrete ocmodel.
%
% CALCFP(OCOBJ) calculates the fix point for the canonical system map of the
% discrete ocmodel OCOBJ using the symbolic toolbox.
%
% CALCFP(OCOBJ,X0) calculates the fix point numerically with starting
% values given at X0. If X0 is empty symbolic calculation is applied. For
% the numerical calculation 'FSOLVE' is used by default.
%
% CALCFP(OCOBJ,X0,ARCARG) a fix point is calculated for the canonical
% systems map corresponding to arc ARCARG. If ARCARG is empty or missing the
% equilibria for every arc are computed.
%
% CALCFP(OCOBJ,X0,ARCARG,ITERATE) a fix point of the ITERATE iterate is calculated.
%
% CALCFP(OCOBJ,X0,ARCARG,OPT) in the ocoption structure OPT options for the
% numerical solver can be provided (category OPT.EQ). In the category
% GENERAL the equation sover can be changed
% opt.GENERAL.EquationSolver='fsolve' or ...
%
% EP = CALCFP(...) the output argument EP is a cell array of 'dynprimitive'
% instances.

depvar0=[]; % initial guess for equilibria
opt=[]; % option for fsolve
numeric=1;
arcarg=[];
iterate=[];

if isempty(docObj)
    return
end
if nargin>=2
    depvar0=varargin{1};
end
if nargin>=3
    arcarg=varargin{2};
end
if nargin>=4
    iterate=varargin{3};
end
if nargin>=5
    opt=varargin{4};
end
if isempty(depvar0)
    % symbolic calculations
    numeric=0;
end
if isempty(opt)
    opt=defaultocoptions;
end
if isempty(iterate)
    iterate=1;
end

if isempty(arcarg)
    arcarg=arcargument(docObj);
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
            numeq=canonicalmapequationnum(docObj,act_arcarg);
            for ii=1:numinitpt
                if numeq*iterate==numcoord
                    epcounter=epcounter+1;
                    [depvar,fval,exitflag,info,jacob]=ocmateqsolve(docObj,'FixPointEquation',[],depvar0(:,ii),act_arcarg,opt,iterate);
                    depvar=reshape(depvar,[],iterate);
                    if iterate>1 && ~all(sum(abs(diff(depvar,[],2)))>1e-6)
                        exitflag=0;
                    end
                    epstruct(epcounter).y=depvar;
                    epstruct(epcounter).x=zeros(1,iterate);
                    epstruct(epcounter).y0=depvar(:,1);
                    epstruct(epcounter).x0=0;
                    epstruct(epcounter).period=iterate;
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
                zerofunc=subsparametervalue(docObj,canonicalmap(docObj,[],act_arcarg,1));
                equationnum=length(zerofunc);
                for ii=1:length(zerofunc)
                    zerofunc(ii)=simplify(zerofunc(ii));
                end
                zerofunc=repmat(zerofunc,iterate,1);
                [dependentvariablenamet dependentvariablenametpn{1:iterate}]=dependentvariable(docObj,act_arcarg,iterate);
                numdependentvar=numel(dependentvariablenamet);
                for ii=2:iterate
                    for jj=1:numdependentvar
                        zerofunc(equationnum*(ii-1)+[1:equationnum])=ocmatsubs(zerofunc(equationnum*(ii-1)+[1:equationnum]),[dependentvariablenametpn{1}{jj} '=' dependentvariablenametpn{ii}{jj}],symkernel);
                        zerofunc(equationnum*(ii-1)+[1:equationnum])=ocmatsubs(zerofunc(equationnum*(ii-1)+[1:equationnum]),[dependentvariablenamet{jj} '=' dependentvariablenametpn{ii-1}{jj}],symkernel);
                    end
                end
                for ii=1:numdependentvar
                    zerofunc(equationnum*(iterate-1)+[1:equationnum])=ocmatsubs(zerofunc(equationnum*(iterate-1)+[1:equationnum]),[dependentvariablenametpn{iterate}{ii} '=' dependentvariablenamet{ii}],symkernel);
                end
                sol=ocmatsolve(prepareequationstring(char(zerofunc),symkernel), ...
                    preparevariablestring(cell2vectorstring([dependentvariablenamet dependentvariablenametpn{1:iterate-1}]),symkernel),symkernel);
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
                            epstruct(epcounter).y(counterdepvar,1)=double(sym(sol(ii).(dependentvariablenamet{counterdepvar})));
                        catch
                            ocmatmsg('Solution value ''%s'' for ''%s'' is not a numeric value.\nSolution discarded.\n',sol(ii).(dependentvariablenamet{counterdepvar}),dependentvariablenamet{counterdepvar})
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
ocFP=cell(epcounter,1);
for ii=1:epcounter
    ocFP{ii}=mapprimitive(epstruct(ii),docObj);
    ocFP{ii}.linearization=jacobian(docObj,ocFP{ii});
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
    case 'mupad'
        varn=['{' removesquarebracket(varn) '}'];
end