function sout=gradcont(varargin)
global OCGRADCONT

[OCGRADCONT.problem_func,sol0,opt]=ParseCommandLine(varargin{:});

delete(findall(0,'Tag','OCMATUserStop')) % delete possible message boxes from a previous continuation

% function handles of the actual continuation type
problemhandles=feval(OCGRADCONT.problem_func);
OCGRADCONT.operatoreq=problemhandles{1};
OCGRADCONT.frechetder=problemhandles{2};
OCGRADCONT.gradientsolution=problemhandles{3};
OCGRADCONT.defaultprocessor=problemhandles{7};
OCGRADCONT.testfunc=problemhandles{8};
OCGRADCONT.targetvaluefunc=problemhandles{9};
OCGRADCONT.probleminit=problemhandles{10};
OCGRADCONT.operatorpfrechet=problemhandles{11};
OCGRADCONT.predictextremal=problemhandles{18};
OCGRADCONT.predictextremaldiff=problemhandles{19};
OCGRADCONT.plotcontinuation=problemhandles{21};
OCGRADCONT.printcontinuation=problemhandles{30};
OCGRADCONT.done=problemhandles{16};
OCGRADCONT.adapt=problemhandles{17};

OCGRADCONT.workspaceadapt=problemhandles{20};
OCGRADCONT.formatsolution=problemhandles{22};
OCGRADCONT.testadmissibility=problemhandles{23};
OCGRADCONT.drearr=problemhandles{24};
OCGRADCONT.saveintermediate=problemhandles{26};
OCGRADCONT.datapath=problemhandles{27};

% define options
OCGRADCONT.OPTIONS.R0=getocoptions(opt,'GRADIENT','PenaltyParameter');
OCGRADCONT.OPTIONS.gradtol=getocoptions(opt,'GRADIENT','GradientTolerance');
OCGRADCONT.OPTIONS.gamma=getocoptions(opt,'GRADIENT','SearchDirectionParameter');
OCGRADCONT.OPTIONS.mu=getocoptions(opt,'GRADIENT','InexactStepSizeParameter');
OCGRADCONT.OPTIONS.maxgraditer=getocoptions(opt,'GRADIENT','MaxGradientIteration');
OCGRADCONT.OPTIONS.maxlinesearchiter=getocoptions(opt,'GRADIENT','MaxLinesearchIteration');
OCGRADCONT.OPTIONS.initlinestepwidth=getocoptions(opt,'GRADIENT','InitLineStepWidth');
OCGRADCONT.OPTIONS.minlinestepwidth=getocoptions(opt,'GRADIENT','MinLineStepWidth');
OCGRADCONT.OPTIONS.maxlinestepwidth=getocoptions(opt,'GRADIENT','MaxLineStepWidth');
OCGRADCONT.OPTIONS.gradientmappinggamma=getocoptions(opt,'GRADIENT','GradientMappingGamma');
OCGRADCONT.OPTIONS.gradientmappingmethod=getocoptions(opt,'GRADIENT','GradientMappingMethod');
OCGRADCONT.OPTIONS.directionalderivativestep=getocoptions(opt,'GRADIENT','DirectionalDerivativeStep');
OCGRADCONT.OPTIONS.projection=getocoptions(opt,'GRADIENT','Projection');
OCGRADCONT.OPTIONS.linestep_increase=getocoptions(opt,'GRADIENT','LineStepIncreaseFactor');
OCGRADCONT.OPTIONS.globalcorrector=getocoptions(opt,'GRADIENT','GlobalCorrector');
OCGRADCONT.OPTIONS.printgradientinfo=getocoptions(opt,'GRADIENT','PrintGradientInfo');

OCGRADCONT.OPTIONS.backward=getocoptions(opt,'OCCONTARG','Backward');
OCGRADCONT.OPTIONS.newtonabstol=getocoptions(opt,'NEWTON','AbsTol');
OCGRADCONT.OPTIONS.newtonreltol=getocoptions(opt,'NEWTON','RelTol');
OCGRADCONT.OPTIONS.maxnewtiter=getocoptions(opt,'NEWTON','MaxNewtonIters');
OCGRADCONT.OPTIONS.maxprobes=getocoptions(opt,'NEWTON','MaxProbes');
OCGRADCONT.OPTIONS.checksingular=logical(getocoptions(opt,'NEWTON','CheckSingular'));
OCGRADCONT.OPTIONS.optimset=opt.EQ;
OCGRADCONT.OPTIONS.maxtestiters=getocoptions(opt,'OCCONTARG','MaxTestIters');
OCGRADCONT.OPTIONS.testtolerance=getocoptions(opt,'OCCONTARG','TestTolerance');
OCGRADCONT.OPTIONS.funtolerance=getocoptions(opt,'OCCONTARG','FunTolerance');
OCGRADCONT.OPTIONS.vartolerance=getocoptions(opt,'OCCONTARG','VarTolerance');
OCGRADCONT.OPTIONS.maxdistance=getocoptions(opt,'OCCONTARG','MaxDistance');

OCGRADCONT.linesearchmethod=getocoptions(opt,'GRADIENT','LineSearchMethod');
OCGRADCONT.gradientsolver=str2func(getocoptions(opt,'GENERAL','GradMethod'));
OCGRADCONT.gradientodesolver=str2func(getocoptions(opt,'GENERAL','ODESolver'));
OCGRADCONT.newtonsolver=str2func(getocoptions(opt,'GENERAL','NewtonSolver'));
OCGRADCONT.OPTIONS.zerodeviationtolerance=getocoptions(opt,'GENERAL','ZeroDeviationTolerance');

OCGRADCONT.OPTIONS.continuationmethod=getocoptions(opt,'OCCONTARG','ContinuationMethod'); % 0 ... Arc-length continuation, 1 Moore-Penrose continuation

MaxContinuationSteps=getocoptions(opt,'OCCONTARG','MaxContinuationSteps');
Backward=getocoptions(opt,'OCCONTARG','Backward');
MakeMovie=getocoptions(opt,'GRADIENT','MakeMovie');
ExitOnTargetValue=getocoptions(opt,'OCCONTARG','ExitOnTargetValue');
PlotCont=logical(strcmp(getocoptions(opt,'OCCONTARG','PlotCont'),'on'));
CheckAdmissibility=getocoptions(opt,'OCCONTARG','CheckAdmissibility');

HitTargetValue=getocoptions(opt,'OCCONTARG','HitTargetValue');
SaveIntermediate=strcmpi(getocoptions(opt,'OCCONTARG','SaveIntermediate'),'on');
WorkSpace=getocoptions(opt,'OCCONTARG','WorkSpace');
PrintContStats=strcmpi(getocoptions(opt,'OCCONTARG','PrintContStats'),'on');

stepwidth=getocoptions(opt,'OCCONTARG','InitStepWidth');
stepwidth_min=getocoptions(opt,'OCCONTARG','MinStepWidth');
stepwidth_max=getocoptions(opt,'OCCONTARG','MaxStepWidth');
stepwidth_inc_fac=getocoptions(opt,'OCCONTARG','IncreaseFactor');
stepwidth_dec_fac=getocoptions(opt,'OCCONTARG','DecreaseFactor');
dir_check_step=getocoptions(opt,'OCCONTARG','CheckStep');
dir_check_angle=getocoptions(opt,'OCCONTARG','CheckAngle');

[coeff0,tangent0,extremal0]=adapt2gradsolver(sol0);


if HitTargetValue
    OCGRADCONT.hitvzero=zeros(2,OCGRADCONT.TargetValueNum);  % tangent where testf is zero
    OCGRADCONT.ahv=1;
end

StartTime = clock;
% initialize user workspace
if WorkSpace
    if OCGRADCONT.probleminit(coeff0,extremal0,tangent0)~=0
        ocmaterror('Initializer failed.');
    end
end

[coeff0,extremal0,tangent0,newtoniter,graditer,dgraditer,lineiter,dlineiter,numrepdd,statinfo]=CorrectStartPoint(coeff0,extremal0,tangent0);
if isempty(coeff0)
    ocmatmsg('No convergence at starting function sol.\n')
    ocmatmsg('elapsed time=%.1f secs\n', etime(clock, StartTime));
    ocmatmsg('0 npoints\n');
    sout=[];
    return;
end
if Backward
    tangent0=-tangent0;
end
s.index=1;
s.label='00';
s.data.sol=OCGRADCONT.formatsolution(coeff0,extremal0,tangent0);
s.msg='This is the first solution of the GRADIENT continuation';

fprintf('first solution found\n');
fprintf('tangent vector to first solution found\n');
fprintf(1,' Mesh size: %g\n',OCGRADCONT.TIMEMESH.num);
fprintf(1,' Newton Iterations: %g\n',newtoniter);
fprintf(1,' Gradient Iterations: %g\n',graditer);
fprintf(1,' Gradient Iterations for Extremal Tangent: %g\n',dgraditer);
fprintf(1,' Extremal Tangent Iterations: %g\n',numrepdd); % maximum number of iterations for successful tangent calculation

sout=s;
if SaveIntermediate
    %extremal0.newtoniterations=newtoniter;
    %extremal0.gradientiterations=graditer;
    contout.extremal=extremal0;
    contout.coeff=coeff0;
    contout.tangent=tangent0;
    if length(extremal0)>1
        contout.objectivevalue=extremal0(1).objectivevalue;
    else
        contout.objectivevalue=extremal0.objectivevalue;
    end
    contout.gradientinfo=statinfo;
    failed=OCGRADCONT.saveintermediate(sout,contout,1);
end
if HitTargetValue
    % WM: calculate all testfunctions at once
    [hitval,failed]=EvalHitFunc(0,coeff0,extremal0,tangent0);
    if ~isempty(hitval)
        OCGRADCONT.hitvzero(2,:)=hitval;
        if ~isempty(failed)
            error('Evaluation of test functions failed at starting solution.');
        end
    else
        HitTargetValue=0;
    end
end

% Open dialog box to interrupt continuation manually
UserStop=stoploop('Stop continuation.');
if PlotCont
    PlotContinuationProcess(coeff0,extremal0,tangent0,contout,MakeMovie);
end
contnum=1;
exitflag=0;
coeff1=coeff0;
extremal1=extremal0;
tangent1=tangent0;
nonadmissible=0;

if MakeMovie
    OCGRADCONT.gifmoviename=[OCGRADCONT.basicmoviefilename '_' OCGRADCONT.problem_func '.gif'];
    actframe=getframe(gcf);
    [imageind,clrmp] = rgb2ind(actframe.cdata,256);
    imwrite(imageind,clrmp,OCGRADCONT.gifmoviename,'gif','DelayTime',0,'Loopcount',inf);
end

while contnum <= MaxContinuationSteps && ~exitflag && ~UserStop.Stop()
    corrections=1;
    while 1 && ~UserStop.Stop()
        [coeffpre,extremalpre]=PredictExtremal(coeff1,extremal1,stepwidth,tangent1);
        [coeff2,extremal2,tangent2,newtoniter,graditer,dgraditer,lineiter,dlineiter,numrepdd,statinfo]=OCGRADCONT.newtonsolver(coeffpre,extremalpre,tangent1);
        if ~isempty(coeff2) && ((contnum < dir_check_step) || ((tangent1'*tangent2 > dir_check_angle))) && lineiter<=OCGRADCONT.OPTIONS.maxlinesearchiter && graditer<=OCGRADCONT.OPTIONS.maxgraditer% && (impr || graditer>1 || lineiter>3) %&& sum((extremalpre.v-extremal2.v).^2)>0% && tangent2(end)>0
            [res,passed]=TestSolutionIntegrity(coeff2,extremal2);
            nonadmissible=0;
            % test admissibility
            if CheckAdmissibility & contnum>=CheckAdmissibility
                [nonadmissible,infoS,labelS]=EvalAdmissibilityFunc(coeff2,extremal2,tangent2);
            end
            if ~nonadmissible
                failed=[];
                if isempty(failed)||~failed
                    break
                end
            end
        end
        if PlotCont
            PlotContinuationProcess(coeff2,extremal2,tangent2,contout,MakeMovie);
        end
        if nonadmissible
            if PrintContStats
                ocmatmsg(' \nNon admissible solution detected, reduce stepwidth.\n')
            end
            s.index=contnum;
            s.label='NAS';
            s.data.sol=OCGRADCONT.formatsolution(coeff1,extremal1,tangent1);
            s.data.soln=OCGRADCONT.formatsolution(coeff2,extremal2,tangent2);
            s.data.infoS=infoS;
            s.data.labelS=labelS;
            s.msg='This is a non-admissible solution.';
            sout=[sout; s];
        end
        if stepwidth > stepwidth_min
            stepwidth=max(stepwidth*stepwidth_dec_fac,stepwidth_min);
            corrections=corrections + 1;
        else      % if not then fail
            ocmatmsg('Current step size too small (point %d)\n',contnum);
            exitflag=1;
            if isempty(coeff2)
                coeff2=coeff1;
                tangent2=tangent1;
                extremal2=extremal1;
            else
                if WorkSpace
                    if OCGRADCONT.probleminit(coeff2,tangent2)~=0
                        ocmaterror('Initializer failed.');
                    end
                end
            end
            break
        end
        if lineiter>OCGRADCONT.OPTIONS.maxlinesearchiter
            if PrintContStats
                ocmatmsg(' \nMaximum number of line iteration steps (%d) exceeded, reduce stepwidth %f.\n',OCGRADCONT.OPTIONS.maxlinesearchiter,stepwidth)
            end
        end
        if graditer>OCGRADCONT.OPTIONS.maxgraditer
            if PrintContStats
                ocmatmsg(' \nMaximum number of gradient steps (%d) exceeded, reduce stepwidth %f.\n',OCGRADCONT.OPTIONS.maxgraditer,stepwidth)
            end
        end
    end
    if PrintContStats
        maxnorm_dv=calcmaxgradientnorm(extremal2);
        fprintf(1,'\n Continuation step No.: %i\n',contnum);
        fprintf(1,' Stepwidth: %g ([%g - %g])\n',stepwidth,stepwidth_min,stepwidth_max);
        fprintf(1,' Mesh size: %g\n',OCGRADCONT.TIMEMESH.num);
        fprintf(1,' Newton Iterations: %g\n',newtoniter);
        fprintf(1,' Maximum Gradient norm: %g\n',maxnorm_dv);
        fprintf(1,' Gradient Iterations: %g\n',graditer);
        fprintf(1,' Line Iterations: %g\n',lineiter);
        fprintf(1,' Difference Prediction-Correction: %g\n',max(abs(coeffpre-coeff2)));
        fprintf(1,' Difference Last new sol: %g\n',max(abs(coeff1-coeff2)));
        fprintf(1,' Gradient Iterations for Extremal Tangent: %g\n',dgraditer);
        fprintf(1,' Extremal Tangent Iterations: %g\n',numrepdd); % maximum number of iterations for successful tangent calculation
        PrintContinuationProcess(coeff2,extremal2,tangent2);
    end
    if stepwidth < stepwidth_max && corrections==1 && newtoniter<5 && graditer<0.3*OCGRADCONT.OPTIONS.maxgraditer %&& numrepdd<=3
        stepwidth=min(stepwidth*stepwidth_inc_fac,stepwidth_max);
    end
    % Hit target value
    if HitTargetValue
        [hitval,failed]=EvalHitFunc(0,coeff2,extremal2,tangent2);
        if isempty(failed) || ~failed
            OCGRADCONT.hitvzero(OCGRADCONT.ahv,:)=hitval;
            OCGRADCONT.ahv=3-OCGRADCONT.ahv;
            testchanges=sign(OCGRADCONT.hitvzero(1,:))~=sign(OCGRADCONT.hitvzero(2,:));
            if any(testchanges)
                testidx=find(testchanges);
                [coeff2,extremal2,tangent2]=FindSolAtTargetPoint(testidx,coeff1,extremal1,tangent1,coeff2,extremal2,tangent2);
                if ~isempty(coeff2)
                    fprintf('\n Target value hit.\n');
                    s.index=contnum;
                    s.label='HTV';
                    s.data.sol=OCGRADCONT.formatsolution(coeff2,extremal2,tangent2);
                    %[failed,s]=feval(OCGRADCONT.process,-1,tmesh2,coeff2,tangent2,s);
                    s.msg='This is the solution at the target value.';
                    sout=[sout; s];
                    %fprintf(' label=%s\n Continuation parameter=%g\n', s.label,coeff2(OCGRADCONT.HE.contparametercoord));
                    if ExitOnTargetValue && contnum>=ExitOnTargetValue
                        exitflag=1;
                    end
                else
                    coeff2=coeff1;
                    tangent2=tangent1;
                    extremal2=extremal1;
                    if ExitOnTargetValue && contnum>=ExitOnTargetValue
                        exitflag=1;
                    end
                end
            end
        end
    end

    if SaveIntermediate
        %extremal2.newtoniterations=newtoniter;
        %extremal2.gradientiterations=graditer;
        contout(contnum+1).extremal=extremal2;
        contout(contnum+1).coeff=coeff2;
        contout(contnum+1).tangent=tangent2;
        if length(extremal2)>1
            contout(contnum+1).objectivevalue=extremal2(1).objectivevalue;
        else
            contout(contnum+1).objectivevalue=extremal2.objectivevalue;
        end
        contout(contnum+1).gradientinfo=statinfo;
        failed=OCGRADCONT.saveintermediate(sout,contout,contnum+1);
    end
    %[b,L]=DiscontinuityDetection(contout,contnum);
%     if PrintContStats
%         fprintf(1,' Continuity Measure: %g\n',L);
%     end
%     if b
%         s.index=contnum;
%         s.label='OVJ';
%         s.data.sol=OCGRADCONT.formatsolution(coeff2,extremal2,tangent2);
%         s.data.L=L;
%         s.msg='Jump in objective value detected.';
%         sout=[sout; s];
%     end
    if PlotCont
        PlotContinuationProcess(coeff2,extremal2,tangent2,contout,MakeMovie);
    end
    if MakeMovie
        actframe=getframe(gcf);
        [imageind,clrmp] = rgb2ind(actframe.cdata,256);
        imwrite(imageind,clrmp,OCGRADCONT.gifmoviename,'gif','DelayTime',0,'WriteMode','append');

    end
    coeff1=coeff2;
    extremal1=extremal2;
    tangent1=tangent2;
    contnum=contnum+1;
    OCGRADCONT.contnum=contnum;
end
fprintf('\n');
EndTime = clock;

s.index=contnum;
s.label='99';
s.data.sol=OCGRADCONT.formatsolution(coeff1,extremal1,tangent1);
s.msg='This is the last solution of the GRADIENT continuation';
sout=[sout; s];

if SaveIntermediate
    failed=OCGRADCONT.saveintermediate(sout,contout,contnum+1);
end

UserStop.Clear();
clear UserStop
fprintf('\n elapsed time  = %.1f secs\n', etime(EndTime, StartTime));


%-----------------------------------
%
% Command line parser
%
%-----------------------------------

function [problemfunc,solinit,opt]=ParseCommandLine(problemfunc,solinit,varargin)

if nargin < 2
    ocmaterror('wrong number of input arguments');
end
if isempty(solinit)
    ocmaterror('provided solution is empty')
end

opt=[];
if nargin > 2
    solinit.solverinfo.tangent=varargin{1};
    if nargin > 3
        opt=varargin{2};
    end
end
if isempty(opt)
    opt=defaultocoptions;
end
%--< END OF CMDL PARSER >--

%----------------------------
%
% Check admissibility of solution
%
%----------------------------

function [failed,infoS,labelS]=EvalAdmissibilityFunc(coeff,extremal,tangent)
global OCGRADCONT

[failed,infoS,labelS]=OCGRADCONT.testadmissibility(coeff,extremal,tangent);

%--< END OF EVALADMISSIBILITYFUNC >--


%---------------------------------
%
% Print the solution
%
%
function b=PrintContinuationProcess(coeff,extremal,tangent)
global OCGRADCONT
b=OCGRADCONT.printcontinuation(coeff,extremal,tangent);

%------------------------------------------
%
%  Evaluate targetvaluefunctions
%
%------------------------------------------

function [out,failed]=EvalHitFunc(id,coeff,extremal,tangent)
global OCGRADCONT

if id==0
    [out,failed]=OCGRADCONT.targetvaluefunc(1:OCGRADCONT.TargetValueNum,coeff,extremal,tangent);
else
    [out,failed]=OCGRADCONT.targetvaluefunc(id,coeff,extremal,tangent);
end

%-------------------------------------
%
% Start point Corrector
%
%-------------------------------------

function [coeff,extremal,tangent,newtoniter,graditer,dgraditer,lineiter,dlineiter,numrepdd,statinfo]=CorrectStartPoint(coeff0,extremal0,tangent0)
global OCGRADCONT

coeff=[];
extremal=[];
tangent=[];
% no tangent vector given, cycle through base-vectors
ndim=length(coeff0);
ii=ndim;
if ~isempty(tangent0)
    [coeff,extremal,tangent,newtoniter,graditer,dgraditer,lineiter,dlineiter,numrepdd,statinfo]=OCGRADCONT.newtonsolver(coeff0,extremal0,tangent0);
end
tangent0=zeros(ndim,1);
while isempty(coeff) && ii>=1
    tangent0(ii)=1;
    try
        [coeff,extremal,tangent,newtoniter,graditer,dgraditer,lineiter,dlineiter,numrepdd,statinfo]=OCGRADCONT.newtonsolver(coeff0,extremal0,tangent0);
    catch
        lasterr
    end
    tangent0(ii)=0;
    ii=ii-1;
end

function [coeffpre,extremalpre]=PredictExtremal(coeff,extremal,stepwidth,tangent)
global OCGRADCONT
[coeffpre,extremalpre]=OCGRADCONT.predictextremal(coeff,extremal,stepwidth,tangent);


%-------------------------------------
%
% Find the solution at the target point
%
%-------------------------------------
function [coeff,extremal,tangent]=FindSolAtTargetPoint(id,coeff1,extremal1,tangent1,coeff2,extremal2,tangent2)
global OCGRADCONT


t1=EvalHitFunc(id,coeff1,extremal1,tangent1);
t2=EvalHitFunc(id,coeff2,extremal2,tangent2);

ii=1;
tangent=[];
coeff=[];
tmax=10*max(abs(t1(id)),abs(t2(id)));
p=1;
while ii<=OCGRADCONT.OPTIONS.maxtestiters
    % WM: make educated guess of where the zero point might be
    if tmax < Inf
        r=abs(t1(id)/(t1(id)-t2(id)))^p;
    else
        r=0.5;
    end
    %coeff3=coeff1+r*(coeff2-coeff1);
    [coeff3,extremal3]=OCGRADCONT.predictextremaldiff(coeff1,coeff2,extremal1,extremal2,r);
    tangent3=tangent1+r*(tangent2-tangent1);
    [coeff,extremal,tangent]=OCGRADCONT.newtonsolver(coeff3,extremal3,tangent3);
    if isempty(coeff)
        coeff=coeff3;
        extremal=extremal3;
        tangent=tangent3;
    end
    tval=EvalHitFunc(id,coeff,extremal,tangent);
    dist1=norm(coeff-coeff1);
    dist2=norm(coeff-coeff2);
    if abs(tval(id))>tmax
        %     fprintf('testfunction behaving badly.\n');
        tangent=[];
        extremal=[];
        coeff=[];
        break;
    end
    if abs(tval(id))<=OCGRADCONT.OPTIONS.testtolerance && min(dist1,dist2) < 100*OCGRADCONT.OPTIONS.vartolerance
        break;
    elseif sign(tval(id))==sign(t2(id))
        coeff2=coeff;
        tangent2=tangent;
        t2(id)=tval(id);
        p=1.02;
    else
        coeff1=coeff;
        tangent1=tangent;
        t1(id)=tval(id);
        p=0.98;
    end
    ii=ii+1;
    coeff=[];
    extremal=[];
    tangent=[];
end


%---------------------------------
%
% Plot the solution
%
%
function b=PlotContinuationProcess(coeff,extremal,tangent,contdata,makemovie)
global OCGRADCONT
b=OCGRADCONT.plotcontinuation(coeff,extremal,tangent,contdata,makemovie);
  


function [res,passed]=TestSolutionIntegrity(coeff,extremal)
global OCGRADCONT
[freepar,modelpar]=OCGRADCONT.drearr(coeff);
passed=1;
res=[];

function [b,L]=DiscontinuityDetection(contout,contnum)
global OCGRADCONT
b=[];
L=NaN;
if contnum<3
    return
end
oval=[contout(end-3:end).objectivevalue];
contpar=[contout(end-3:end).coeff];
contpar=contpar(end,:);

ointerp=interp1(contpar(1:3),oval(1:3),contpar(end),'cubic','extrap');
L=abs((ointerp-oval(end))/oval(end));
b=L>OCGRADCONT.OPTIONS.objectivejumpthreshold;

function nrm=calcmaxgradientnorm(extremal)
global OCGRADCONT
nrm=-inf;

for ii=1:length(extremal)
    if ~isfield(extremal(ii),'dv')
        return
    end
    nrm=max(nrm,sqrt(max(sum((extremal(ii).dv.^2)))));
end