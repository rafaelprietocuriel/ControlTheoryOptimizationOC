function sout=cont_static(varargin)
global OCSTATCONT

[OCSTATCONT.problem_func,sol0,opt]=ParseCommandLine(varargin{:});

delete(findall(0,'Tag','OCMATUserStop')) % delete possible message boxes from a previous continuation

% Open dialog box to interrupt continuation manually
UserStop=stoploop('Stop continuation.');

% function handles of the actual continuation type
problemhandles=feval(OCSTATCONT.problem_func);
OCSTATCONT.generalizedoperatoreq=problemhandles{1};
OCSTATCONT.generalizedfrechetder=problemhandles{2};
OCSTATCONT.defaultprocessor=problemhandles{7};
OCSTATCONT.testfunc=problemhandles{8};
OCSTATCONT.targetvaluefunc=problemhandles{9};
OCSTATCONT.probleminit=problemhandles{10};
OCSTATCONT.plotcontinuation=problemhandles{21};
OCSTATCONT.printcontinuation=problemhandles{30};
OCSTATCONT.singmat=problemhandles{13};
OCSTATCONT.process=problemhandles{14};
OCSTATCONT.locate=problemhandles{15};
OCSTATCONT.done=problemhandles{16};
OCSTATCONT.adapt=problemhandles{17};
try
    OCSTATCONT.dataadaptation=problemhandles{18};
catch
    OCSTATCONT.dataadaptation=[];
end
OCSTATCONT.workspaceadapt=problemhandles{20};
OCSTATCONT.formatsolution=problemhandles{22};
OCSTATCONT.drearr=problemhandles{24};
OCSTATCONT.saveintermediate=problemhandles{26};

% define options
OCSTATCONT.OPTIONS.generalizedderivativethreshold=getocoptions(opt,'NEWTON','GeneralizedDerivativeThreshold');
OCSTATCONT.OPTIONS.newtonabstol=getocoptions(opt,'NEWTON','AbsTol');
OCSTATCONT.OPTIONS.newtonreltol=getocoptions(opt,'NEWTON','RelTol');
OCSTATCONT.OPTIONS.maxnewtiter=getocoptions(opt,'NEWTON','MaxNewtonIters');
OCSTATCONT.OPTIONS.maxprobes=getocoptions(opt,'NEWTON','MaxProbes');
OCSTATCONT.OPTIONS.checksingular=logical(getocoptions(opt,'NEWTON','CheckSingular'));
OCSTATCONT.OPTIONS.singularthreshold=getocoptions(opt,'NEWTON','SingularThreshold');
OCSTATCONT.OPTIONS.maxtestiters=getocoptions(opt,'OCCONTARG','MaxTestIters');
OCSTATCONT.OPTIONS.testtolerance=getocoptions(opt,'OCCONTARG','TestTolerance');
OCSTATCONT.OPTIONS.funtolerance=getocoptions(opt,'OCCONTARG','FunTolerance');
OCSTATCONT.OPTIONS.vartolerance=getocoptions(opt,'OCCONTARG','VarTolerance');

OCSTATCONT.newtonsolver=str2func(getocoptions(opt,'GENERAL','NewtonSolver'));

OCSTATCONT.OPTIONS.continuationmethod=getocoptions(opt,'OCCONTARG','ContinuationMethod'); % 0 ... Arc-length continuation, 1 Moore-Penrose continuation

Backward=getocoptions(opt,'OCCONTARG','Backward');
stepwidth=getocoptions(opt,'OCCONTARG','InitStepWidth');
stepwidth_min=getocoptions(opt,'OCCONTARG','MinStepWidth');
stepwidth_max=getocoptions(opt,'OCCONTARG','MaxStepWidth');
stepwidth_inc_fac=getocoptions(opt,'OCCONTARG','IncreaseFactor');
stepwidth_dec_fac=getocoptions(opt,'OCCONTARG','DecreaseFactor');
dir_check_step=getocoptions(opt,'OCCONTARG','CheckStep');
dir_check_angle=getocoptions(opt,'OCCONTARG','CheckAngle');

MaxContinuationSteps=getocoptions(opt,'OCCONTARG','MaxContinuationSteps');
MakeMovie=getocoptions(opt,'GRADIENT','MakeMovie');
ExitOnTargetValue=getocoptions(opt,'OCCONTARG','ExitOnTargetValue');
PlotCont=logical(strcmp(getocoptions(opt,'OCCONTARG','PlotCont'),'on'));
CheckAdmissibility=getocoptions(opt,'OCCONTARG','CheckAdmissibility');

HitTargetValue=getocoptions(opt,'OCCONTARG','HitTargetValue');
SaveIntermediate=strcmpi(getocoptions(opt,'OCCONTARG','SaveIntermediate'),'on');
WorkSpace=getocoptions(opt,'OCCONTARG','WorkSpace');
PrintContStats=strcmpi(getocoptions(opt,'OCCONTARG','PrintContStats'),'on');

[coeff0,tangent0]=adapt2statsolver(sol0);


if HitTargetValue
    OCSTATCONT.hitvzero=zeros(2,OCSTATCONT.TargetValueNum);  % tangent where testf is zero
    OCSTATCONT.ahv=1;
end


StartTime = clock;
% initialize user workspace
if WorkSpace
    if OCSTATCONT.probleminit(coeff0,tangent0)~=0
        ocmaterror('Initializer failed.');
    end
end

[coeff0,tangent0,newtoniter]=CorrectStartPoint(coeff0,tangent0);
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
s.data.sol=OCSTATCONT.formatsolution(coeff0,tangent0);
s.msg='This is the first solution of the OPTIM continuation';

fprintf(1,'first solution found\n');
fprintf(1,'tangent vector to first solution found\n');
fprintf(1,' Newton Iterations: %g\n',newtoniter);

sout=s;
if SaveIntermediate
    contout.coeff=coeff0;
    contout.tangent=tangent0;
    failed=OCSTATCONT.saveintermediate(sout,contout,1);
end
if HitTargetValue
    % WM: calculate all testfunctions at once
    [hitval,failed]=EvalHitFunc(0,coeff0,tangent0);
    if ~isempty(hitval)
        OCSTATCONT.hitvzero(2,:)=hitval;
        if ~isempty(failed)
            error('Evaluation of test functions failed at starting solution.');
        end
    else
        HitTargetValue=0;
    end
end
contnum=1;
exitflag=0;
coeff1=coeff0;
tangent1=tangent0;

if MakeMovie
    OCSTATCONT.gifmoviename=[OCSTATCONT.basicmoviefilename '_' OCSTATCONT.problem_func '.gif'];
    actframe=getframe(gcf);
    [imageind,clrmp] = rgb2ind(actframe.cdata,256);
    imwrite(imageind,clrmp,OCSTATCONT.gifmoviename,'gif','DelayTime',0,'Loopcount',inf);
end

while contnum <= MaxContinuationSteps && ~exitflag && ~UserStop.Stop()
    corrections=1;
    while 1 && ~UserStop.Stop()
        coeffpre=coeff1+stepwidth*tangent1;
        [coeff2,tangent2,newtoniter]=OCSTATCONT.newtonsolver(coeffpre,tangent1);
        if ~isempty(coeff2) && ((contnum < dir_check_step) || ((tangent1'*tangent2 > dir_check_angle))) % && (impr || graditer>1 || lineiter>3) %&& sum((extremalpre.v-extremal2.v).^2)>0% && tangent2(end)>0
            [res,passed]=TestSolutionIntegrity(coeff2);
            nonadmissible=0;
            % test admissibility
            if CheckAdmissibility & contnum>=CheckAdmissibility
                %[nonadmissible,infoS,labelS]=EvalAdmissibilityFunc(coeff2,extremal2,tangent2);
            end
            if ~nonadmissible
                failed=[];
                if isempty(failed)||~failed
                    break
                end
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
            else
                if WorkSpace
                    if OCSTATCONT.probleminit(coeff2,tangent2)~=0
                        ocmaterror('Initializer failed.');
                    end
                end
            end
            break
        end
    end
    if PrintContStats
        fprintf(1,'\n Continuation step No.: %i\n',contnum);
        fprintf(1,' stepwidth: %g ([%g - %g])\n',stepwidth,stepwidth_min,stepwidth_max);
        fprintf(1,' Newton Iterations: %g\n',newtoniter);
        PrintContinuationProcess(coeff2,tangent2);
    end
    if stepwidth < stepwidth_max && newtoniter<5
        stepwidth=min(stepwidth*stepwidth_inc_fac,stepwidth_max);
    end
    if HitTargetValue
        [hitval,failed]=EvalHitFunc(0,coeff2,tangent2);
        if isempty(failed) || ~failed
            OCSTATCONT.hitvzero(OCSTATCONT.ahv,:)=hitval;
            OCSTATCONT.ahv=3-OCSTATCONT.ahv;
            testchanges=sign(OCSTATCONT.hitvzero(1,:))~=sign(OCSTATCONT.hitvzero(2,:));
            if any(testchanges)
                testidx=find(testchanges);
                [coeff2,tangent2]=FindSolAtTargetPoint(testidx,coeff1,tangent1,coeff2,tangent2);
                if ~isempty(coeff2)
                    fprintf('\n Target value hit.\n');
                    s.index=contnum;
                    s.label='HTV';
                    s.data.sol=OCSTATCONT.formatsolution(coeff2,tangent2);
                    %[failed,s]=feval(OCSTATCONT.process,-1,tmesh2,coeff2,tangent2,s);
                    s.msg='This is the solution at the target value.';
                    sout=[sout; s];
                    %fprintf(' label=%s\n Continuation parameter=%g\n', s.label,coeff2(OCSTATCONT.HE.contparametercoord));
                    if ExitOnTargetValue && contnum>=ExitOnTargetValue
                        exitflag=1;
                    end
                else
                    coeff2=coeff1;
                    tangent2=tangent1;
                    if ExitOnTargetValue && contnum>=ExitOnTargetValue
                        exitflag=1;
                    end
                end
            end
        end
    end
    if SaveIntermediate
        contout(contnum+1).coeff=coeff1;
        contout(contnum+1).tangent=tangent1;
        failed=OCSTATCONT.saveintermediate(sout,contout,contnum+1);
    end
    if PlotCont
        PlotContinuationProcess(coeff2,tangent2,contout,MakeMovie);
    end
    if MakeMovie
        actframe=getframe(gcf);
        [imageind,clrmp] = rgb2ind(actframe.cdata,256);
        imwrite(imageind,clrmp,OCSTATCONT.gifmoviename,'gif','DelayTime',0,'WriteMode','append');

    end
    coeff1=coeff2;
    tangent1=tangent2;
    contnum=contnum+1;
    OCSTATCONT.contnum=contnum;
end
fprintf('\n');
EndTime = clock;

s.index=contnum;
s.label='99';
s.data.sol=OCSTATCONT.formatsolution(coeff1,tangent1);
s.msg='This is the last solution of the GRADIENT continuation';
sout=[sout; s];

if SaveIntermediate
    failed=OCSTATCONT.saveintermediate(sout,contout,contnum+1);
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
function b=PrintContinuationProcess(coeff,tangent)
global OCSTATCONT
b=OCSTATCONT.printcontinuation(coeff,tangent);

%------------------------------------------
%
%  Evaluate targetvaluefunctions
%
%------------------------------------------

function [out,failed]=EvalHitFunc(id,coeff,tangent)
global OCSTATCONT

if id==0
    [out,failed]=OCSTATCONT.targetvaluefunc(1:OCSTATCONT.TargetValueNum,coeff,tangent);
else
    [out,failed]=OCSTATCONT.targetvaluefunc(id,coeff,tangent);
end

%-------------------------------------
%
% Start point Corrector
%
%-------------------------------------

function [coeff,tangent,newtoniter]=CorrectStartPoint(coeff0,tangent0)
global OCSTATCONT
coeff=[];
tangent=[];
% no tangent vector given, cycle through base-vectors
ndim=length(coeff0);
ii=ndim;
if ~isempty(tangent0)
    [coeff,tangent,newtoniter]=OCSTATCONT.newtonsolver(coeff0,tangent0);
end
tangent0=zeros(ndim,1);
while isempty(coeff) && ii>=1
    tangent0(ii)=1;
    try
        [coeff,tangent,newtoniter]=OCSTATCONT.newtonsolver(coeff0,tangent0);
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
function [coeff,tangent]=FindSolAtTargetPoint(id,coeff1,tangent1,coeff2,tangent2)
global OCSTATCONT


t1=EvalHitFunc(id,coeff1,tangent1);
t2=EvalHitFunc(id,coeff2,tangent2);

ii=1;
tangent=[];
coeff=[];
tmax=10*max(abs(t1(id)),abs(t2(id)));
p=1;
while ii<=OCSTATCONT.OPTIONS.maxtestiters
    % WM: make educated guess of where the zero point might be
    if tmax < Inf
        r=abs(t1(id)/(t1(id)-t2(id)))^p;
    else
        r=0.5;
    end
    coeff3=coeff1+r*(coeff2-coeff1);
    tangent3=tangent1+r*(tangent2-tangent1);
    [coeff,tangent]=OCSTATCONT.newtonsolver(coeff3,tangent3);
    if isempty(coeff)
        coeff=coeff3;
        tangent=tangent3;
    end
    tval=EvalHitFunc(id,coeff,tangent);
    dist1=norm(coeff-coeff1);
    dist2=norm(coeff-coeff2);
    if abs(tval(id))>tmax
        %     fprintf('testfunction behaving badly.\n');
        tangent=[];
        coeff=[];
        break;
    end
    if abs(tval(id))<=OCSTATCONT.OPTIONS.testtolerance && min(dist1,dist2) < 100*OCSTATCONT.OPTIONS.vartolerance
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
    tangent=[];
end


%---------------------------------
%
% Plot the solution
%
%
function b=PlotContinuationProcess(coeff,tangent,contdata,makemovie)
global OCSTATCONT
b=OCSTATCONT.plotcontinuation(coeff,tangent,contdata,makemovie);


function [res,passed]=TestSolutionIntegrity(coeff,extremal)
global OCSTATCONT
[freepar,modelpar]=OCSTATCONT.drearr(coeff);
passed=1;
res=[];

