function varargout=ppdecont(ppdeObj,InitBVStruct,varargin)
%
% OCCONT continue an extremal solution of an oc model
%
% OCXP=OCCONT(OCOBJ,INITSTRUCT) continuation of a solution
% specified in INITSTRUCT of the oc model OCOBJ
% OCXP is the last regularly computed ocextremal of the continuation
% proccess. INITSTRUCT is a structur returned by INITOCCONT(...).
%
% OCXP=OCCONT(CONTARG,OCOBJ,INITSTRUCT,OPT) OPT provides ocoptions for the
% continuation process

% global variable
if exist('OcBVPVar','var')
    clear global OcBVPVar
end
global OcBVPVar
delete(findall(0,'Tag','OCMATUserStop')) % delete possible message boxes from a previous continuation

OcBVPVar=InitBVStruct;
occontargument(-1);
occontargument(InitBVStruct);
bvsolactual=[];
bvsolnonadmissible=[];
intermediateresult=[];
relativedistance=[];
stepwidthhistory=[];
objectivevalue=[];
if ~isstruct(InitBVStruct) || ~isfield(InitBVStruct,'BVPFunction') || isempty(InitBVStruct.BVPFunction)
    return
end

%ocXP=eval([class(initbvsol) '([])']);

if isempty(ppdeObj)
    return
end

opt=[];
conttype='';
if nargin>=3
    opt=varargin{1};
end

if isempty(opt)
    opt=defaultocoptions;
end
% initialize options
jacfuncflag=strcmp(getocoptions(opt,'OCCONTARG','BVPFJacobian'),'on');
bvpjacfuncflag=strcmp(getocoptions(opt,'OCCONTARG','BVPBCJacobian'),'on');
bvpmaxnum=getocoptions(opt,'OCCONTARG','BVPMaxNum');
bvpreltol=getocoptions(opt,'BVP','RelTol');
bvpabstol=getocoptions(opt,'BVP','AbsTol');
vectorizeflag=getocoptions(opt,'BVP','Vectorized');
bvpsolver=getocoptions(opt,'GENERAL','BVPMethod');
if isfield(OcBVPVar.InitBVSolution,'solver') && ~strcmp(OcBVPVar.InitBVSolution.solver,bvpsolver)
    OcBVPVar.InitBVSolution=rmfield(OcBVPVar.InitBVSolution,'solver');
end
initbvsol=OcBVPVar.InitBVSolution;

maxcontnum=getocoptions(opt,'OCCONTARG','MaxContinuationSteps');
maxgridpoints=getocoptions(opt,'OCCONTARG','MaxGridPoints');
initstepwidth=getocoptions(opt,'OCCONTARG','InitStepWidth');
contstepwidth=getocoptions(opt,'OCCONTARG','StepWidth');
maxresnorm=getocoptions(opt,'OCCONTARG','MaxResiduumNorm');
totreldist=getocoptions(opt,'OCCONTARG','TotalRelativeDistance');
meaniter=getocoptions(opt,'OCCONTARG','MeanIteration');
minstep=getocoptions(opt,'OCCONTARG','MinStepWidth');
maxstep=getocoptions(opt,'OCCONTARG','MaxStepWidth');
backward=getocoptions(opt,'OCCONTARG','Backward');
checkflag=strcmp(getocoptions(opt,'OCCONTARG','Check'),'on');
checktimeflag=strcmp(getocoptions(opt,'OCCONTARG','CheckTime'),'on');
checkstart=getocoptions(opt,'OCCONTARG','CheckStart');
admTol=getocoptions(opt,'OC','AdmissibleTolerance');
incfact=getocoptions(opt,'OCCONTARG','StepFactor');

occontargument('FSolveOptions',getfield(opt,'EP'));
OcBVPVar.FSolveOptions=getfield(opt,'EP');
tst=occontargument('BVPFunction.Continuation');
% Argument passed to the BVP file
OcBVPVar.ActContinuationType=0; % 0 ... fixed stepwidth width, 1 ... adaptive quadratic, 2 ... adaptive linear
OcBVPVar.PreviousBVSol=[];
OcBVPVar.StepWidth=contstepwidth;


% initbvsol.parameters=0;
% OcBVPVar.PathFollowing.Index=2;
% OcBVPVar.PathFollowing.maxcont=50;
% OcBVPVar.ActContinuationType=-10;
% OcBVPVar.PathFollowing.Tangent=initbvsol.tangent;

% initialize variables
conttype=InitBVStruct.ContinuationType;
continuationarg=InitBVStruct.ContinuationArgument;
contid=lower(InitBVStruct.Identifier);
outdir=fullfile(getocmatpath('data'),ppdeObj.ocInfos.Name);
if backward
    OcBVPVar.ContinuationVector=-OcBVPVar.ContinuationVector;
end
jacfunc='';
bvpjacfunc='';
ppfunc='';
blockincreasenum=0;
bvsolnonadmissible=[];

% initialize display and save options
displayflag=strcmp(getocoptions(opt,'OCCONTARG','Display'),'on');
saveflag=strcmp(getocoptions(opt,'OCCONTARG','SaveOutput'),'on');
saveintermediateflag=strcmp(getocoptions(opt,'OCCONTARG','SaveIntermediate'),'on');
if displayflag % generate axis labels
    switch contid
        case {'xtendtime','i2eps','i2epws','i2epsts','i2epu','i2epss','i2epuu','i2epo','boundary', ...
                'i2epsc','i2epcs','i2eps_i2eps','i2epcs_i2eps','i2epsc_i2eps', ...
                'i2eps_i2epu','i2eps_i2epuu','fin','fini','fini_fini','i2lc', ...
                'ev','freendtime','initpec','criticalmanifold','xtendtimecriticalmanifold', ...
                'freendtime_freendtime','h2eps','userdynamics','i2dhts','i2dhts_i2dhts','i2epwu'}
            if 1%~transformid(1)
                %OcBVPVar.StateLabel=eval(['{' var2symvar(makestr('x_',1:ppdeObj.ocProblem.State.Num+ppdeObj.ocProblem.UncontrolledState.Num)) '}']);
                OcBVPVar.StateLabel=eval(['{' var2symvar(makestr('x_',1:ppdeObj.ocProblem.State.Num)) '}']);
                OcBVPVar.DynVarLabel=eval(['{' var2symvar(makestr('\lambda_',1:ppdeObj.ocProblem.State.Num)) '}']);
                OcBVPVar.ContParLabel={'contpar'};
            end
        case {'p2ep','p2eps','pclc'}
            if 1%~transformid(1)
                OcBVPVar.ContParLabel={'contpar'};
                OcBVPVar.LimSetLabel=eval(['{' var2symvar(makestr('E_',1:ppdeObj.ocProblem.State.Num)) '}']);
                OcBVPVar.StateLabel=eval(['{' var2symvar(makestr('x_',1:ppdeObj.ocProblem.State.Num)) '}']);
                OcBVPVar.DynVarLabel=eval(['{' var2symvar(makestr('\lambda_',1:ppdeObj.ocProblem.State.Num)) '}']);
            end
        otherwise
            warning('''%s'' not implemented yet. Display is disabled.',contid)
            displayflag=0;
    end
    if strcmp(getocoptions(opt,'OCCONTARG','DisplayWindow'),'docked')
        set(gcf,'WindowStyle','docked')
    end
end

actualreldistance=0;

bcfunc=str2func(OcBVPVar.BVPFunction.General);
odefunc=str2func(OcBVPVar.ODEFunction.General);
if jacfuncflag
    jacfunc=str2func(OcBVPVar.ODEFunction.GeneralJacobian);
end
if bvpjacfuncflag
    bvpjacfunc=str2func(['Jac4' lower(InitBVStruct.Identifier)]);
end
if isfield(OcBVPVar.BVPFunction,'PostProcessing')
    ppfunc=str2func(OcBVPVar.BVPFunction.PostProcessing);
end

%bvpoptions=bvpset('BCJacobian', bvpjacfunc,'FJacobian',jacfunc,'RelTol', bvpreltol, 'AbsTol',bvpabstol, 'NMax', bvpmaxnum, 'Stats', 'off','Vectorized',vectorizeflag);
%bvpoptions=opt.BVP;

contflag=1; % do the continuation while contflag
initialpoints=[];
initialpointsnonadmissible=[];
limitpoints=[];
endpoints=[];
icontswitch=0;

contnum=0; % counts the number of continuation steps
stepwidth=initstepwidth;
laststep=0; % indicating if the last step is reached
bvsol=cell(3,1);
bvsol{3}=initbvsol; % bvsol{3} denotes the actual solution

%actualreldistance=stepwidth; % set the distance for the fixed step problem
%OcBVPVar.InitialValue=OcBVPVar.StartValue+actualreldistance*OcBVPVar.ContinuationVector;
% Open dialog box to interrupt continuation manually
UserStop=stoploop('Stop continuation.');

while contflag && ~UserStop.Stop()
    fprintf(1, '\n Continuation step No.: %i\n',contnum);
    fprintf(1, ' stepwidth: %g\n',stepwidth);

    % store solutions from the two previous continuation steps
    failnum=0;
    while 1
        % initialize continuation type
        if ~laststep && strcmp(conttype,'aq') && contnum==2 && ~error_flag
            stepwidth=contstepwidth;
            OcBVPVar.ActContinuationType=1;
            if ~isfield(bvsol{3},'parameters')
                bvsol{2}.parameters=actualreldistance-stepwidth;
                bvsol{3}.parameters=actualreldistance;
            else
                bvsol{2}.parameters(end+1)=actualreldistance-stepwidth;
                bvsol{3}.parameters(end+1)=actualreldistance;
            end
        end
        if ~laststep && strcmp(conttype,'al') && contnum==2 && ~error_flag
            stepwidth=contstepwidth;
            OcBVPVar.ActContinuationType=2;
            if ~isfield(bvsol{3},'parameters')
                bvsol{2}.parameters=actualreldistance-stepwidth;
                bvsol{3}.parameters=actualreldistance;
            else
                bvsol{2}.parameters(end+1)=actualreldistance-stepwidth;
                bvsol{3}.parameters(end+1)=actualreldistance;
            end
        end
        if contnum > 1
            % approximating the continued solution from the two previous
            % steps
            try
                bvsolpredict=[];
                bvsolpredict.x=bvsol{3}.x;
                try
                    bvsolpredict.y = bvsol{3}.y*(1+alphafac) - alphafac*deval(bvsol{2},bvsol{3}.x);  % predicted solution
                catch
                    bvsolpredict.y = bvsol{3}.y*(1+alphafac) - alphafac*deval6c(bvsol{2},bvsol{3}.x);  % predicted solution
                end
                if isfield(bvsol{3},'parameters') && isfield(bvsol{2},'parameters')
                    if ~laststep
                        bvsolpredict.parameters = bvsol{3}.parameters(:).'*(1+alphafac) - alphafac*bvsol{2}.parameters(:).';
                    else
                        bvsolpredict.parameters = bvsol{3}.parameters(:).'*(1+alphafac) - alphafac*bvsol{2}.parameters(1:length(bvsol{3}.parameters)).';
                    end
                end
                if OcBVPVar.ActContinuationType>=2
                    OcBVPVar.PreviousBVSol=[[bvsol{2}.y(1:OcBVPVar.CanSysDim,1); bvsol{2}.parameters(end)],[bvsol{3}.y(1:OcBVPVar.CanSysDim,1); bvsol{3}.parameters(end)]];
                elseif OcBVPVar.ActContinuationType==1
                    OcBVPVar.PreviousBVSol=[bvsol{3}.y(:,1); bvsol{3}.parameters(:)];
                end
                %bvsolpredict=bvsol{3};
            catch
                bvsolpredict=bvsol{3};
            end
        else
            bvsolpredict = bvsol{3};
            stepwidth=initstepwidth;
        end
        alphafac=1;
        try
            solver
            if isfield(OcBVPVar,'NewModelParameterValue')
                ppdeObj.ocInstant.ParameterValue=OcBVPVar.NewModelParameterValue;
            end
            if ~error_flag && checkflag &&  contnum>=checkstart
                if strcmpi(OcBVPVar.CheckingType,'unique')
                    bvsoltest=bvsolactual;
                    if isfield(OcBVPVar,'Coordinate')
                        if iscell(OcBVPVar.Coordinate)
                            bvsoltest.y=bvsolactual.y([OcBVPVar.Coordinate{:}],:);
                        else
                            bvsoltest.y=bvsolactual.y(OcBVPVar.Coordinate,:);
                        end
                    end
                    if strcmp(OcBVPVar.CalcObjectiveValue,'on')
                        bvsoltest.y([OcBVPVar.ObjectiveValueCoordinate{:}],:)=[];
                    end

                    if strcmp(OcBVPVar.ContinuationArgument,{'heteroclinic','homoclinic'})
                        par=OcBVPVar.ModelParameterValue;
                        par(OcBVPVar.TargetCoordinate)=OcBVPVar.InitialValue;
                        par(OcBVPVar.FreeParameter)=bvsoltest.parameters(end-length(OcBVPVar.FreeParameter)+1:end);
                        testocObj=changeparameter(ppdeObj,par);
                    else
                        testocObj=ppdeObj;
                        testocObj.ocInstant.ParameterValue=OcBVPVar.ModelParameterValue;
                    end
                    [violationflag violationmat violationterminal,userviolation]=testadmissibility4unique(testocObj,bvsoltest,OcBVPVar.ArcIdentifier,OcBVPVar.Function.Admissible,opt,OcBVPVar.TerminalInequalityConstraint,OcBVPVar.Function.UserConstraint,OcBVPVar.ObjectiveDynamicsCoordinates,OcBVPVar.ExcludeConstraint);
                    if ~isempty(OcBVPVar.Function.UserConstraint)
                        violationflag=violationflag||any(userviolation);
                    end
                elseif  strcmpi(OcBVPVar.CheckingType,'multiple')
                    [violationflag violationmat]=testadmissibility4multiple(ppdeObj,bvsolactual,OcBVPVar.ArcIdentifier,opt);
                end

                if checktimeflag
                    % check positivity of timeintervals
                    [tviolationflag tviolationmat]=testtimeintervals(bvsoltest,OcBVPVar,opt);
                else
                    tviolationflag=0;

                end

                if violationflag || tviolationflag
                    warning('Calculated trajectory is not admissible.')
                    bvsolnonadmissible=bvsolactual;
                    bvsolnonadmissible.violation=violationmat;
                    error_flag=1;
                end
            else
                bvsoltest=bvsolactual;
                % check positivity of timeintervals
                %[tviolationflag tviolationmat]=testtimeintervals(bvsoltest,OcBVPVar,opt);
                tviolationflag=0;
                if tviolationflag
                    warning('Timeinterval(s) of calculated trajectory is not admissible.')
                    bvsolnonadmissible=bvsolactual;
                    bvsolnonadmissible.violation=violationmat;
                    error_flag=1;
                end

            end
        catch
            lasterr
            error_flag=1;
        end

        if ~error_flag % solution found
            temporary2constant; % replace values by temporarily calculated values during the bvp process
            if OcBVPVar.ActContinuationType
                actualreldistance=bvsolactual.parameters(end);
            end
            fprintf(1, ' Relative distance: %g\n',totreldist-actualreldistance)
            if iternum<=meaniter && contnum>2 && blockincreasenum>=0 && contnum-icontswitch-2>0 % increase stepwidth width
                if incfact*stepwidth<maxstep
                    stepwidth=stepwidth*incfact;
                    alphafac = incfact*alphafac;
                    icontswitch=contnum;
                else
                    alphafac=maxstep/stepwidth;
                    stepwidth = maxstep;
                end
            end
            if blockincreasenum<0
                blockincreasenum=blockincreasenum+1;
            end

            if contnum>maxcontnum  || laststep || stepwidth<minstep
                contflag=0;
            end

            if strcmp(OcBVPVar.CalcObjectiveValue,'on')
                objectivevalue=[objectivevalue bvsolactual.y(([OcBVPVar.ObjectiveValueCoordinate{:}]),end)];
            end
            if ~OcBVPVar.ActContinuationType
                if totreldist-actualreldistance-stepwidth<0
                    alphafac=(totreldist-actualreldistance)/stepwidth;
                    stepwidth=totreldist-actualreldistance;
                    laststep=1;
                    if isfield(bvsolactual,'parameters')
                        initialpoints=[initialpoints [bvsolactual.y(:,1);bvsolactual.parameters(:);actualreldistance]];
                        endpoints=[endpoints [bvsolactual.y(:,end);bvsolactual.parameters(:);actualreldistance]];
                    else
                        initialpoints=[initialpoints [bvsolactual.y(:,1);actualreldistance]];
                        endpoints=[endpoints [bvsolactual.y(:,end);actualreldistance]];
                    end

                    if isfield(OcBVPVar,'SaddlePoint')
                        if ~iscell(OcBVPVar.SaddlePoint)
                            limitpoints=[limitpoints OcBVPVar.SaddlePoint];
                        else
                            limitpoints=[limitpoints vertcat(OcBVPVar.SaddlePoint{:})];
                        end
                    end
                end
            else
                if totreldist-actualreldistance<0
                    stepwidth=totreldist-bvsol{3}.parameters(end);
                    actualreldistance=bvsol{3}.parameters(end);
                    alphafac=1;%stepwidth/(bvsol{3}.parameters(end)-bvsol{2}.parameters(end));
                    laststep=1;
                    OcBVPVar.ActContinuationType=0;
                    if length(bvsol{3}.parameters)>1
                        bvsol{3}.parameters(end)=[];
                    else
                        bvsol{3}=rmfield(bvsol{3},'parameters');
                    end
                    bvsolactual=bvsol{3};
                    %bvsol=bvsol([1 2 2]);
                end
            end

            if ~laststep
                if OcBVPVar.ActContinuationType%isfield(bvsolactual,'parameters')
                    %initialpoints=[initialpoints [bvsolactual.y(1:OcBVPVar.ArcDimension,1);OcBVPVar.ArcIdentifier(1);transformid(1)]];
                    initialpoints=[initialpoints [bvsolactual.y(:,1);bvsolactual.parameters(:)]];
                    endpoints=[endpoints [bvsolactual.y(:,end);bvsolactual.parameters(:)]];
                else
                    %initialpoints=[initialpoints [bvsolactual.y(1:OcBVPVar.ArcDimension,1);OcBVPVar.ArcIdentifier(1);transformid(1)]];
                    if isfield(bvsolactual,'parameters')
                        initialpoints=[initialpoints [bvsolactual.y(:,1);bvsolactual.parameters(:);actualreldistance]];
                        endpoints=[endpoints [bvsolactual.y(:,end);bvsolactual.parameters(:);actualreldistance]];
                    else
                        initialpoints=[initialpoints [bvsolactual.y(:,1);actualreldistance]];
                        endpoints=[endpoints [bvsolactual.y(:,end);actualreldistance]];
                    end
                end
                if isfield(OcBVPVar,'SaddlePoint')
                    if ~iscell(OcBVPVar.SaddlePoint)
                        limitpoints=[limitpoints OcBVPVar.SaddlePoint];
                    else
                        limitpoints=[limitpoints vertcat(OcBVPVar.SaddlePoint{:})];
                    end
                end
            end
            % display and store solution
            if displayflag
                plotbvsol(initialpoints,limitpoints,bvsolactual,endpoints)
            end
            if ~OcBVPVar.ActContinuationType
                actualreldistance=actualreldistance+stepwidth; % set the distance for the fixed step problem
                %OcBVPVar.InitialValue=actualreldistance;
                OcBVPVar.InitialValue=OcBVPVar.StartValue+actualreldistance*OcBVPVar.ContinuationVector;
            else
                OcBVPVar.StepWidth=stepwidth;
            end

            % update cell of solutions
            bvsol(1)=[];
            bvsol{3}=bvsolactual;

            break;         % pursue next continuation step
        else
            if laststep
                laststep=0;
            end
            if failnum>5 % interrupt continuation after five succesive failures
                warning(['Continuation process interrupted since the number of successive failures exceeds five.'])
                contflag=0;
                break
            end
            failnum=failnum+1;
            if ~OcBVPVar.ActContinuationType
                if ~isempty(bvsolnonadmissible)
                    if isfield(bvsolactual,'parameters')
                        %initialpointsnonadmissible=[initialpointsnonadmissible [bvsolnonadmissible.y(:,1);bvsolnonadmissible.parameters(:);actualreldistance]];
                    else
                        %initialpointsnonadmissible=[initialpointsnonadmissible [bvsolnonadmissible.y(:,1);bvsolnonadmissible]];
                    end
                end
                actualreldistance=actualreldistance-stepwidth; %  reduce distance to its previous value
            end
            blockincreasenum=-25;
            %blockincreasenum=-5;
            stepwidth=stepwidth/incfact;
            alphafac = alphafac/incfact;
            if stepwidth<minstep
                warning(['Continuation process interrupted since minimum stepwidth is reached.'])
                contflag=0;
                break
            end
            if ~OcBVPVar.ActContinuationType
                actualreldistance=actualreldistance+stepwidth; % set the distance for the fixed step problem
                %OcBVPVar.InitialValue=actualreldistance;
                OcBVPVar.InitialValue=OcBVPVar.StartValue+actualreldistance*OcBVPVar.ContinuationVector;
            else
                OcBVPVar.StepWidth=stepwidth;
            end
            continue;
        end
    end

    if saveflag
        savebvsol
    end
    contnum=contnum+1;
end
UserStop.Clear();
clear UserStop

if nargout>=1
    if ~isempty(bvsolactual)
        if isocasymptotic(OcBVPVar.InitOcTrajectory)
            if isfield(OcBVPVar,'SlowFastSwitch')
                if isfield(OcBVPVar,'InitialTimes') && length(OcBVPVar.InitialTimes)==1
                    inputarg={'slowfastswitch',OcBVPVar.SlowFastSwitch};
                else
                    inputarg={'slowfastswitch',OcBVPVar.SlowFastSwitch,'timehorizon',OcBVPVar.TimeHorizon};
                end
            else
                inputarg{1}=[];
            end
            varargout{1}=ode2trj(OcBVPVar.InitOcTrajectory,bvsolactual,[],OcBVPVar.Jacobian,OcBVPVar.SaddlePoint,'timehorizon',OcBVPVar.TimeHorizon,inputarg{:});
        else
            varargout{1}=ode2trj(OcBVPVar.InitOcTrajectory,bvsolactual,[],'timehorizon',OcBVPVar.TimeHorizon);
        end
    else
        varargout{1}=[];
    end
end
if nargout>=2
    if ~isempty(bvsolnonadmissible)
        if isocasymptotic(OcBVPVar.InitOcTrajectory)
            varargout{2}=ode2trj(OcBVPVar.InitOcTrajectory,bvsolnonadmissible,bvsolnonadmissible.violation,OcBVPVar.Jacobian,OcBVPVar.SaddlePoint);
        else
            varargout{2}=ode2trj(OcBVPVar.InitOcTrajectory,bvsolnonadmissible,bvsolnonadmissible.violation,'timehorizon',OcBVPVar.TimeHorizon);
        end
    else
        varargout{2}=[];
    end
end

if nargout>=3
    varargout{3}=bvsolactual;

end
if nargout>=4
    if strncmpi(OcBVPVar.Identifier,'p',1)
        warning('Parameter value changes. Actual ocmodel is replaced. Stored information is lost.')
        varargout{3}=changeparameter(ppdeObj,OcBVPVar.TargetCoordinate,OcBVPVar.InitialValue);
    else
        varargout{4}=ppdeObj;
    end
end

% Sub functions
    function solver
        iternum=0;
        lastwarn('');
        switch bvpsolver
            case 'bvp4c'
                bvsolactual= bvp4c(odefunc,bcfunc,bvsolpredict,opt.BVP);
                % if error is thrown for missing fields 'bcevals'
                % check if the original m-file is adapted to return the
                % number of iterations 'bcevals' see the documentation of
                % OCMat
            case 'bvp5c'
                bvsolactual= bvp5c(odefunc,bcfunc,bvsolpredict,opt.BVP);
            case 'bvp6c'
                bvsolactual= bvp6c(odefunc,bcfunc,bvsolpredict,opt.BVP);
            case 'bvp6cdefault'
                bvsolactual= bvp6cdefault(odefunc,bcfunc,bvsolpredict,opt.BVP);
            case 'sbvp'
                bvsolactual=bvsolpredict;
                [bvsolactual.x,bvsolactual.y,bvsolactual.xcol,bvsolactual.ycol,bcres]= sbvp(str2func(OcBVPVar.SBVPFunction.General),bvsolpredict.x,bvsolpredict.y,opt.SBVP);
                bvsolactual.bcevals=10;
            otherwise
                error('BV solve ''%s'' not implemented.',bvpsolver)
        end
        if ~isempty(ppfunc)
            feval(ppfunc,bvsolactual);
        end
        if ~OcBVPVar.ActContinuationType
            bvsolactual.contparameter=actualreldistance;
        else
            bvsolactual.contparameter=bvsolactual.parameters(end);
        end
        switch bvpsolver
            case 'sbvp'

            otherwise
                % returns the residuum of the actual solution
                if isfield(bvsolactual,'parameters')
                    bcres=feval(bcfunc, bvsolactual.y(:,1), bvsolactual.y(:,end),bvsolactual.parameters);
                else
                    bcres=feval(bcfunc, bvsolactual.y(:,1), bvsolactual.y(:,end));
                end
        end

        fprintf(1, ' Mesh size: %d\n', size(bvsolactual.y,2));
        if isfield(bvsolactual,'bcevals')
            fprintf(1, ' Iteration number: %d\n',bvsolactual.bcevals);
            iternum=bvsolactual.bcevals;
        else
            warning('The number of calls to the BC function is not returned. (See the ''readme.1st'' file).')
            iternum=inf;
        end
        fprintf(1, ' Residuum norm: %g\n',norm(bcres(1:end)))
        if norm(bcres(1:end)) > 100*maxresnorm || length(bvsolactual.x)>maxgridpoints || ~isempty(lastwarn)
            error_flag = 1;
        else
            error_flag = 0;
        end
    end

    function savebvsol
        if saveintermediateflag
            intermediateresult=[intermediateresult {bvsolactual}];
            %if ~OcBVPVar.ActContinuationType
            relativedistance=[relativedistance actualreldistance];
            stepwidthhistory=[stepwidthhistory OcBVPVar.StepWidth];
            %end
        end
        save(fullfile(outdir,'ContinuationData'),'bvsolactual','initialpoints','InitBVStruct','bvsolnonadmissible','OcBVPVar','intermediateresult','relativedistance','stepwidthhistory','initialpointsnonadmissible','objectivevalue')
    end
end

function temporary2constant()

global OcBVPVar

if ~isfield(OcBVPVar,'Temporary')
    return
end

fn=fieldnames(OcBVPVar.Temporary);
for ii=1:length(fn)
    if ~isempty(OcBVPVar.Temporary.(fn{ii}))
        OcBVPVar.(fn{ii})=OcBVPVar.Temporary.(fn{ii});
        %OcBVPVar.Temporary.(fn{ii})=[];
    end
end
end