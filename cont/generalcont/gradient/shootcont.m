function sout=shootcont(varargin)
global OCSHOOTCONT

[OCSHOOTCONT.problem_func,sol0,opt]=ParseCommandLine(varargin{:});

delete(findall(0,'Tag','OCMATUserStop')) % delete possible message boxes from a previous continuation

% function handles of the actual continuation type
problemhandles=feval(OCSHOOTCONT.problem_func);
OCSHOOTCONT.operatoreq=problemhandles{1};
OCSHOOTCONT.frechetder=problemhandles{2};
OCSHOOTCONT.gradientsolution=problemhandles{3};
OCSHOOTCONT.defaultprocessor=problemhandles{7};
OCSHOOTCONT.testfunc=problemhandles{8};
OCSHOOTCONT.targetvaluefunc=problemhandles{9};
OCSHOOTCONT.probleminit=problemhandles{10};
OCSHOOTCONT.operatorpfrechet=problemhandles{11};
OCSHOOTCONT.predictextremal=problemhandles{18};
OCSHOOTCONT.predictextremaldiff=problemhandles{19};
OCSHOOTCONT.plotcontinuation=problemhandles{21};
OCSHOOTCONT.printcontinuation=problemhandles{30};
OCSHOOTCONT.done=problemhandles{16};
OCSHOOTCONT.adapt=problemhandles{17};

OCSHOOTCONT.workspaceadapt=problemhandles{20};
OCSHOOTCONT.formatsolution=problemhandles{22};
OCSHOOTCONT.drearr=problemhandles{24};
OCSHOOTCONT.saveintermediate=problemhandles{26};
OCSHOOTCONT.datapath=problemhandles{27};

% define options

OCSHOOTCONT.OPTIONS.backward=getocoptions(opt,'OCCONTARG','Backward');
OCSHOOTCONT.OPTIONS.newtonabstol=getocoptions(opt,'NEWTON','AbsTol');
OCSHOOTCONT.OPTIONS.newtonreltol=getocoptions(opt,'NEWTON','RelTol');
OCSHOOTCONT.OPTIONS.maxnewtiter=getocoptions(opt,'NEWTON','MaxNewtonIters');
OCSHOOTCONT.OPTIONS.maxprobes=getocoptions(opt,'NEWTON','MaxProbes');
OCSHOOTCONT.OPTIONS.checksingular=logical(getocoptions(opt,'NEWTON','CheckSingular'));
OCSHOOTCONT.OPTIONS.optimset=opt.EQ;
OCSHOOTCONT.OPTIONS.maxtestiters=getocoptions(opt,'OCCONTARG','MaxTestIters');
OCSHOOTCONT.OPTIONS.testtolerance=getocoptions(opt,'OCCONTARG','TestTolerance');
OCSHOOTCONT.OPTIONS.funtolerance=getocoptions(opt,'OCCONTARG','FunTolerance');
OCSHOOTCONT.OPTIONS.vartolerance=getocoptions(opt,'OCCONTARG','VarTolerance');

OCSHOOTCONT.shootingsolver=str2func(getocoptions(opt,'GENERAL','ShootMethod'));
OCSHOOTCONT.shootingodesolver=str2func(getocoptions(opt,'GRADIENT','ODESolver'));
OCSHOOTCONT.newtonsolver=str2func(getocoptions(opt,'GENERAL','NewtonSolver'));

OCSHOOTCONT.OPTIONS.continuationmethod=getocoptions(opt,'OCCONTARG','ContinuationMethod'); % 0 ... Arc-length continuation, 1 Moore-Penrose continuation
OCSHOOTCONT.OPTIONS.correctcontroltype='prior';
OCSHOOTCONT.OPTIONS.directionalderivativestep=1e-5;

MaxContinuationSteps=getocoptions(opt,'OCCONTARG','MaxContinuationSteps');
stepwidth_min=getocoptions(opt,'OCCONTARG','MinStepWidth');
Backward=getocoptions(opt,'OCCONTARG','Backward');
MakeMovie=getocoptions(opt,'GRADIENT','MakeMovie');
ExitOnTargetValue=getocoptions(opt,'OCCONTARG','ExitOnTargetValue');
PlotCont=logical(strcmp(getocoptions(opt,'OCCONTARG','PlotCont'),'on'));

HitTargetValue=getocoptions(opt,'OCCONTARG','HitTargetValue');
SaveIntermediate=strcmpi(getocoptions(opt,'OCCONTARG','SaveIntermediate'),'on');
WorkSpace=getocoptions(opt,'OCCONTARG','WorkSpace');
PrintContStats=strcmpi(getocoptions(opt,'OCCONTARG','PrintContStats'),'on');

stepwidth=getocoptions(opt,'OCCONTARG','InitStepWidth');
stepwidth_max=getocoptions(opt,'OCCONTARG','MaxStepWidth');
stepwidth_inc_fac=getocoptions(opt,'OCCONTARG','IncreaseFactor');
stepwidth_dec_fac=getocoptions(opt,'OCCONTARG','DecreaseFactor');
dir_check_step=getocoptions(opt,'OCCONTARG','CheckStep');
dir_check_angle=getocoptions(opt,'OCCONTARG','CheckAngle');

[coeff0,tangent0,extremal0]=adapt2shootsolver(sol0);


if HitTargetValue
    OCSHOOTCONT.hitvzero=zeros(2,OCSHOOTCONT.TargetValueNum);  % tangent where testf is zero
    OCSHOOTCONT.ahv=1;
end

StartTime = clock;
% initialize user workspace
if WorkSpace
    if OCSHOOTCONT.probleminit(coeff0,extremal0,tangent0)~=0
        ocmaterror('Initializer failed.');
    end
end

[coeff0,extremal0,tangent0,newtoniter]=CorrectStartPoint(coeff0,extremal0,tangent0);
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
s.data.sol=OCSHOOTCONT.formatsolution(coeff0,extremal0,tangent0);
s.msg='This is the first solution of the GRADIENT continuation';

fprintf('first solution found\n');
fprintf('tangent vector to first solution found\n');
fprintf(1,' Mesh size: %g\n',OCSHOOTCONT.TIMEMESH.num);
fprintf(1,' Newton Iterations: %g\n',newtoniter);

sout=s;
if SaveIntermediate
    contout.extremal=extremal0;
    contout.coeff=coeff0;
    contout.tangent=tangent0;
    failed=OCSHOOTCONT.saveintermediate(sout,contout,1);
end
if HitTargetValue
    % WM: calculate all testfunctions at once
    [hitval,failed]=EvalHitFunc(0,coeff0,extremal0,tangent0);
    if ~isempty(hitval)
        OCSHOOTCONT.hitvzero(2,:)=hitval;
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
    delete(gcf)
    PlotContinuationProcess(coeff0,extremal0,tangent0,contout,MakeMovie);
end
contnum=1;
exitflag=0;
coeff1=coeff0;
extremal1=extremal0;
tangent1=tangent0;

if MakeMovie
    OCSHOOTCONT.gifmoviename=[OCSHOOTCONT.basicmoviefilename '_' OCSHOOTCONT.problem_func '.gif'];
    actframe=getframe(gcf);
    [imageind,clrmp] = rgb2ind(actframe.cdata,256);
    imwrite(imageind,clrmp,OCSHOOTCONT.gifmoviename,'gif','DelayTime',0,'Loopcount',inf);
end

while contnum <= MaxContinuationSteps && ~exitflag && ~UserStop.Stop()
    corrections=1;
    while 1 && ~UserStop.Stop()
        coeffpre=coeff1+stepwidth*tangent1;
        [coeff2,extremal2,tangent2,newtoniter]=OCSHOOTCONT.newtonsolver(coeffpre,extremal1,tangent1);
        if ~isempty(coeff2) && ((contnum < dir_check_step) || ((tangent1'*tangent2 > dir_check_angle)))
            [res,passed]=TestSolutionIntegrity(coeff2,extremal2);
            if passed
                failed=[];
            end
            if isempty(failed)||~failed
                break
            end
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
                    if OCSHOOTCONT.probleminit(coeff2,tangent2)~=0
                        ocmaterror('Initializer failed.');
                    end
                end
            end
            break
        end
    end
    if PrintContStats
        fprintf(1,'\n Continuation step No.: %i\n',contnum);
        fprintf(1,' Stepwidth: %g ([%g - %g])\n',stepwidth,stepwidth_min,stepwidth_max);
        fprintf(1,' Mesh size: %g\n',OCSHOOTCONT.TIMEMESH.num);
        fprintf(1,' Newton Iterations: %g\n',newtoniter);
        PrintContinuationProcess(coeff2,extremal2,tangent2);
    end
    if stepwidth < stepwidth_max && corrections==1 && newtoniter<5 %&& numrepdd<=3
        stepwidth=min(stepwidth*stepwidth_inc_fac,stepwidth_max);
    end
    % Hit target value
    if HitTargetValue
        [hitval,failed]=EvalHitFunc(0,coeff2,extremal2,tangent2);
        if isempty(failed) || ~failed
            OCSHOOTCONT.hitvzero(OCSHOOTCONT.ahv,:)=hitval;
            OCSHOOTCONT.ahv=3-OCSHOOTCONT.ahv;
            testchanges=sign(OCSHOOTCONT.hitvzero(1,:))~=sign(OCSHOOTCONT.hitvzero(2,:));
            if any(testchanges)
                testidx=find(testchanges);
                [coeff2,extremal2,tangent2]=FindSolAtTargetPoint(testidx,coeff1,extremal1,tangent1,coeff2,extremal2,tangent2);
                if ~isempty(coeff2)
                    fprintf('\n Target value hit.\n');
                    s.index=contnum;
                    s.label='HTV';
                    s.data.sol=OCSHOOTCONT.formatsolution(coeff2,extremal2,tangent2);
                    %[failed,s]=feval(OCSHOOTCONT.process,-1,tmesh2,coeff2,tangent2,s);
                    s.msg='This is the solution at the target value.';
                    sout=[sout; s];
                    %fprintf(' label=%s\n Continuation parameter=%g\n', s.label,coeff2(OCSHOOTCONT.HE.contparametercoord));
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
        contout(contnum+1).extremal=extremal2;
        contout(contnum+1).coeff=coeff2;
        contout(contnum+1).tangent=tangent2;
        contout(contnum+1).stepwidth=stepwidth;
        failed=OCSHOOTCONT.saveintermediate(sout,contout,contnum+1);
    end
    if PlotCont
        PlotContinuationProcess(coeff2,extremal2,tangent2,contout,MakeMovie);
    end
    if MakeMovie
        actframe=getframe(gcf);
        [imageind,clrmp] = rgb2ind(actframe.cdata,256);
        imwrite(imageind,clrmp,OCSHOOTCONT.gifmoviename,'gif','DelayTime',0,'WriteMode','append');

    end
    coeff1=coeff2;
    extremal1=extremal2;
    tangent1=tangent2;
    contnum=contnum+1;
    OCSHOOTCONT.contnum=contnum;
end
fprintf('\n');
EndTime = clock;

s.index=contnum;
s.label='99';
s.data.sol=OCSHOOTCONT.formatsolution(coeff1,extremal1,tangent1);
s.msg='This is the last solution of the GRADIENT continuation';
sout=[sout; s];

if SaveIntermediate
    failed=OCSHOOTCONT.saveintermediate(sout,contout,contnum+1);
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

%---------------------------------
%
% Print the solution
%
%
function b=PrintContinuationProcess(coeff,extremal,tangent)
global OCSHOOTCONT
b=OCSHOOTCONT.printcontinuation(coeff,extremal,tangent);

%------------------------------------------
%
%  Evaluate targetvaluefunctions
%
%------------------------------------------

function [out,failed]=EvalHitFunc(id,coeff,extremal,tangent)
global OCSHOOTCONT

if id==0
    [out,failed]=OCSHOOTCONT.targetvaluefunc(1:OCSHOOTCONT.TargetValueNum,coeff,extremal,tangent);
else
    [out,failed]=OCSHOOTCONT.targetvaluefunc(id,coeff,extremal,tangent);
end

%-------------------------------------
%
% Start point Corrector
%
%-------------------------------------

function [coeff,extremal,tangent,newtoniter]=CorrectStartPoint(coeff0,extremal0,tangent0)
global OCSHOOTCONT

coeff=[];
extremal=[];
tangent=[];
% no tangent vector given, cycle through base-vectors
ndim=length(coeff0);
ii=ndim;
if ~isempty(tangent0)
    [coeff,extremal,tangent,newtoniter]=OCSHOOTCONT.newtonsolver(coeff0,extremal0,tangent0);
end
tangent0=zeros(ndim,1);
while isempty(coeff) && ii>=1
    tangent0(ii)=1;
    try
        [coeff,extremal,tangent,newtoniter]=OCSHOOTCONT.newtonsolver(coeff0,extremal0,tangent0);
    catch
        lasterr
    end
    tangent0(ii)=0;
    ii=ii-1;
end


%-------------------------------------
%
% Find the solution at the target point
%
%-------------------------------------
function [coeff,extremal,tangent]=FindSolAtTargetPoint(id,coeff1,extremal1,tangent1,coeff2,extremal2,tangent2)
global OCSHOOTCONT


t1=EvalHitFunc(id,coeff1,extremal1,tangent1);
t2=EvalHitFunc(id,coeff2,extremal2,tangent2);

ii=1;
tangent=[];
coeff=[];
tmax=10*max(abs(t1(id)),abs(t2(id)));
p=1;
while ii<=OCSHOOTCONT.OPTIONS.maxtestiters
    % WM: make educated guess of where the zero point might be
    if tmax < Inf
        r=abs(t1(id)/(t1(id)-t2(id)))^p;
    else
        r=0.5;
    end
    %coeff3=coeff1+r*(coeff2-coeff1);
    [coeff3,extremal3]=OCSHOOTCONT.predictextremaldiff(coeff1,coeff2,extremal1,extremal2,r);
    tangent3=tangent1+r*(tangent2-tangent1);
    [coeff,extremal,tangent]=OCSHOOTCONT.newtonsolver(coeff3,extremal3,tangent3);
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
    if abs(tval(id))<=OCSHOOTCONT.OPTIONS.testtolerance && min(dist1,dist2) < 100*OCSHOOTCONT.OPTIONS.vartolerance
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
global OCSHOOTCONT
b=OCSHOOTCONT.plotcontinuation(coeff,extremal,tangent,contdata,makemovie);
  


function [res,passed]=TestSolutionIntegrity(coeff,extremal)
global OCSHOOTCONT
[freepar,modelpar]=OCSHOOTCONT.drearr(coeff);
passed=1;
res=[];
