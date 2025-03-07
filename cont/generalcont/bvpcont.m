function sout=bvpcont(varargin)
%
% BVPCONT main OCMat continuation file for BVPs
%
% BVPCONT(PROBLEMTYPE,INITSOL) argument PROBLEMTYPE is a string
% characterizing the type of problem to be solved:
%   'extremal2ep'   ... saddle-path of an equilibrium, continuing along the
%                       initial point
%   'extremalp2ep'  ... saddle-path of an equilibrium, continuing along a
%                       parameter value
%   'extremalt2ep'  ... saddle-path of an equilibrium, continuing the
%                       truncation time
%   'indifferencesolution'  ... continuation of an indifference threshold
%                       (Skiba point)
%   'limitextremal' ... continuation of a limitpoint solution
%
% INITSOL is an initial function structure to start the continuation,
% returned by an initialization function, e.g. initocmat_AE_EP,
% initocmat_AE_IS (see OCMat manual)
%
% BVPCONT(PROBLEMTYPE,INITSOL,INITTANGENT) INITTANGENT is an initial
% tangent, usually this argument is empty.
%
% BVPCONT(PROBLEMTYPE,INITSOL,INITTANGENT,OPT) the option structure OPT provides a
% multitude of settings for the continuation process and BVP solver (see
% OCMat manual).
%
% SOUT=BVPCONT(...) the structure array consists at least of two elements
% from the intial and last step of the continuation.
%
% The main structure of the bvpcont file is taken from the MatCont file 'cont'

global OCMATCONT
[OCMATCONT.problem_func,sol0,opt]=ParseCommandLine(varargin{:});
delete(findall(0,'Tag','OCMATUserStop')) % delete possible message boxes from a previous continuation

% function handles of the actual continuation type
problemhandles=feval(OCMATCONT.problem_func);
OCMATCONT.operatoreq=problemhandles{1};
OCMATCONT.frechetder=problemhandles{2};
OCMATCONT.options=problemhandles{3};
OCMATCONT.defaultprocessor=problemhandles{7};
OCMATCONT.testfunc=problemhandles{8};
OCMATCONT.targetvaluefunc=problemhandles{9};
OCMATCONT.probleminit=problemhandles{10};
OCMATCONT.operatorpfrechet=problemhandles{11};
OCMATCONT.plotcontinuation=problemhandles{21};
OCMATCONT.printcontinuation=problemhandles{30};
OCMATCONT.singmat=problemhandles{13};
OCMATCONT.process=problemhandles{14};
OCMATCONT.locate=problemhandles{15};
OCMATCONT.done=problemhandles{16};
OCMATCONT.adapt=problemhandles{17};
try
    OCMATCONT.dataadaptation=problemhandles{18};
catch
    OCMATCONT.dataadaptation=[];
end
OCMATCONT.stateadapt=problemhandles{19};
OCMATCONT.workspaceadapt=problemhandles{20};
OCMATCONT.formatsolution=problemhandles{22};
OCMATCONT.testadmissibility=problemhandles{23};
OCMATCONT.drearr=problemhandles{24};
OCMATCONT.saveintermediate=problemhandles{26};
OCMATCONT.domaindiscretization=problemhandles{28};

ode=problemhandles{4}{1};
bc=problemhandles{4}{2};
if getocoptions(opt,'BVP','FJacobian')
    odejac=problemhandles{5}{1};
else
    odejac=[];
end
if getocoptions(opt,'BVP','BCJacobian')
    bcjac=problemhandles{5}{2};
else
    bcjac=[];
end
if 0%getocoptions(opt,'SBVPOC','FDerivative2')
    try
        odehess=problemhandles{5}{3};
    catch
        odehess=[];
    end
else
    odehess=[];
end
if 0%getocoptions(opt,'SBVPOC','BCDerivative2')
    try
        bchess=problemhandles{5}{4};
    catch
        bchess=[];
    end
else
    bchess=[];
end
if 0%getocoptions(opt,'SBVPOC','FDerivative3')
    try
        odetensor3=problemhandles{5}{5};
    catch
        odetensor3=[];
    end
else
    odetensor3=[];
end
if 0%getocoptions(opt,'SBVPOC','BCDerivative3')
    try
        bctensor3=problemhandles{5}{6};
    catch
        bctensor3=[];
    end
else
    bctensor3=[];
end
icfunjac=[];
icfun=[];


OCMATCONT.OPTIONS.newtonabstol=getocoptions(opt,'NEWTON','AbsTol');
OCMATCONT.OPTIONS.newtonreltol=getocoptions(opt,'NEWTON','RelTol');
OCMATCONT.OPTIONS.maxnewtiter=getocoptions(opt,'NEWTON','MaxNewtonIters');
OCMATCONT.OPTIONS.maxprobes=getocoptions(opt,'NEWTON','MaxProbes');
OCMATCONT.OPTIONS.trm=getocoptions(opt,'NEWTON','TRM');
OCMATCONT.OPTIONS.lambdamin=getocoptions(opt,'NEWTON','LambdaMin');
OCMATCONT.OPTIONS.updatejacfactor=getocoptions(opt,'NEWTON','UpdateJacFactor');
OCMATCONT.OPTIONS.switchtoffnfactor=getocoptions(opt,'NEWTON','SwitchToFFNFactor');
OCMATCONT.OPTIONS.checksingular=logical(getocoptions(opt,'NEWTON','CheckSingular'));
OCMATCONT.OPTIONS.display=getocoptions(opt,'NEWTON','Display');
OCMATCONT.OPTIONS.log=getocoptions(opt,'NEWTON','Log');
OCMATCONT.OPTIONS.singularthreshold=getocoptions(opt,'NEWTON','SingularThreshold');
OCMATCONT.newtonsolver=str2func(getocoptions(opt,'GENERAL','NewtonSolver'));
OCMATCONT.OPTIONS.admissibletol=getocoptions(opt,'GENERAL','AdmissibleTolerance');
OCMATCONT.OPTIONS.zerotimedifftolerance=getocoptions(opt,'GENERAL','ZeroTimeDiffTolerance');

%OCMATCONT.OPTIONS.residualreductionguard=getocoptions(opt,'SBVPOC','residualReductionGuard');  % to prevent mesh oscillations.
OCMATCONT.OPTIONS.meshadaptabstol=getocoptions(opt,'SBVPOC','MeshAdaptAbsTol');
OCMATCONT.OPTIONS.meshadaptreltol=getocoptions(opt,'SBVPOC','MeshAdaptRelTol');
%OCMATCONT.OPTIONS.meshadaptabstol=getocoptions(opt,'NEWTON','AbsTol');
%OCMATCONT.OPTIONS.meshadaptreltol=getocoptions(opt,'NEWTON','RelTol');
OCMATCONT.OPTIONS.xyvectorized=strcmpi(getocoptions(opt,'BVP','Vectorized'),'on');
%OCMATCONT.OPTIONS.meshadaptk=getocoptions(opt,'SBVPOC','MeshAdaptK');
%OCMATCONT.OPTIONS.finemesh=getocoptions(opt,'SBVPOC','MeshAdaptFineMesh');
%OCMATCONT.OPTIONS.meshadaptmaxiter=getocoptions(opt,'SBVPOC','MeshAdaptMaxIter');
OCMATCONT.OPTIONS.maxgridnum=getocoptions(opt,'BVP','Nmax');
%OCMATCONT.OPTIONS.maxnewpoints=getocoptions(opt,'SBVPOC','MaxNewPts');

OCMATCONT.OPTIONS.totalrelativedistance=getocoptions(opt,'OCCONTARG','TotalRelativeDistance');
OCMATCONT.OPTIONS.maxtestiters=getocoptions(opt,'OCCONTARG','MaxTestIters');
OCMATCONT.OPTIONS.testtolerance=getocoptions(opt,'OCCONTARG','TestTolerance');
OCMATCONT.OPTIONS.funtolerance=getocoptions(opt,'OCCONTARG','FunTolerance');
OCMATCONT.OPTIONS.increment=getocoptions(opt,'OCCONTARG','Increment');
OCMATCONT.OPTIONS.vartolerance=getocoptions(opt,'OCCONTARG','VarTolerance');
OCMATCONT.OPTIONS.saveintermediate=strcmp(getocoptions(opt,'OCCONTARG','SaveIntermediate'),'on');
OCMATCONT.OPTIONS.continuationmethod=getocoptions(opt,'OCCONTARG','ContinuationMethod'); % 0 ... Arc-length continuation, 1 Moore-Penrose continuation
OCMATCONT.OPTIONS.maxdistance=getocoptions(opt,'OCCONTARG','MaxDistance');
OCMATCONT.OPTIONS.maxdistancecoordinate=getocoptions(opt,'OCCONTARG','MaxDistanceCoordinate');

OCMATCONT.OPTIONS.messagedisplay=strcmpi(getocoptions(opt,'INIT','MessageDisplay'),'on');
OCMATCONT.OPTIONS.testmodus=strcmpi(getocoptions(opt,'INIT','TestModus'),'on');

OCMATCONT.bvpmethod=getocoptions(opt,'GENERAL','BVPMethod');
OCMATCONT.OPTIONS.multipointbvp=false;
% if ~strcmp(OCMATCONT.bvpmethod,'bvp5c')
%     OCMATCONT.OPTIONS.multipointbvp=logical(strcmp(getocoptions(opt,'GENERAL','MultiPointBVP'),'on'));
% else
%     OCMATCONT.OPTIONS.multipointbvp=false;
% end

% Handle argument functions and additional arguments
[tmesh0,coeff0,tangent0]=adapt2solver(sol0,ode,odejac,bc,bcjac,icfun,icfunjac,odehess,bchess,odetensor3,bctensor3);

HitTargetValue=getocoptions(opt,'OCCONTARG','HitTargetValue');
Singularities=getocoptions(opt,'OCCONTARG','Singularities');
SaveIntermediate=strcmpi(getocoptions(opt,'OCCONTARG','SaveIntermediate'),'on');
PrintContStats=strcmpi(getocoptions(opt,'OCCONTARG','PrintContStats'),'on');
PlotCont=strcmpi(getocoptions(opt,'OCCONTARG','PlotCont'),'on');
WorkSpace=getocoptions(opt,'OCCONTARG','WorkSpace');
IgnoreSings=getocoptions(opt,'OCCONTARG','IgnoreSingularity');
ExitOnTargetValue=getocoptions(opt,'OCCONTARG','ExitOnTargetValue');
Locators=getocoptions(opt,'OCCONTARG','Locators');
Userfunctions=getocoptions(opt,'OCCONTARG','Userfunctions');
MaxContStep=getocoptions(opt,'OCCONTARG','MaxContinuationSteps');
Backward=getocoptions(opt,'OCCONTARG','Backward');
CheckAdmissibility=getocoptions(opt,'OCCONTARG','CheckAdmissibility');
Adapt=getocoptions(opt,'OCCONTARG','Adapt');
MeshAdaptation=getocoptions(opt,'SBVPOC','MeshAdaptation');

% used options in the continuation file
stepwidth=getocoptions(opt,'OCCONTARG','InitStepWidth'); % initialize step width
stepwidth_max=getocoptions(opt,'OCCONTARG','MaxStepWidth');
stepwidth_min=getocoptions(opt,'OCCONTARG','MinStepWidth');
stepwidth_inc_fac=getocoptions(opt,'OCCONTARG','IncreaseFactor');
stepwidth_dec_fac=getocoptions(opt,'OCCONTARG','DecreaseFactor');
dir_check_step=getocoptions(opt,'OCCONTARG','CheckStep');
dir_check_angle=getocoptions(opt,'OCCONTARG','CheckAngle');

OCMATCONT.nActTest=0;
if Singularities
    [OCMATCONT.S,SingLables]=feval(OCMATCONT.singmat);
    [OCMATCONT.nSing, OCMATCONT.nTest]=size(OCMATCONT.S);
    Singularities=OCMATCONT.nSing>0 & OCMATCONT.nTest>0;
    % setup testfunction variables and stuff
    %

    if Singularities
        % Ignore Singularities
        OCMATCONT.S(IgnoreSings,:)=8;
        ActSing=setdiff(1:OCMATCONT.nSing, IgnoreSings);
        OCMATCONT.ActSing=ActSing;
        OCMATCONT.nActSing=length(ActSing);

        % Which test functions must vanish somewhere in S?
        ActTest=find( sum((OCMATCONT.S==0),1) > 0 );
        OCMATCONT.ActTest=ActTest;
        OCMATCONT.nActTest=length(ActTest);
        % WM: Build matrix with indices of testfunctions which
        % should be zero for each singularity
        OCMATCONT.SZ=zeros(OCMATCONT.nTest+1,OCMATCONT.nActSing+1);
        ml=2;
        for si=ActSing(1:OCMATCONT.nActSing)
            t=find( OCMATCONT.S(si,:)==0 )';
            l=size(t,1);
            OCMATCONT.SZ(1:l,si)=t;
            ml=max(ml,l);
        end
        OCMATCONT.SZ=OCMATCONT.SZ(1:ml,:);
        OCMATCONT.atv=1;
        OCMATCONT.testvals=zeros(2, OCMATCONT.nActTest);
        % 1st row: testvals at x1, 2nd: testvals at x2
        OCMATCONT.testzerosol(OCMATCONT.nActTest).tmesh=[];
        OCMATCONT.testzerosol(OCMATCONT.nActTest).coeff=[];
        OCMATCONT.testzerosol(OCMATCONT.nActTest).tangent=[];
        if isempty(ActTest)
            Singularities=0;
        end
    end
end

if HitTargetValue
    OCMATCONT.hitvzero=zeros(2,OCMATCONT.TargetValueNum);  % tangent where testf is zero
    OCMATCONT.ahv=1;
end

if PlotCont
    clf
end
StartTime = clock;

% initialize user workspace
if WorkSpace
    if OCMATCONT.probleminit(tmesh0,coeff0,tangent0)~=0
        ocmaterror('Initializer failed.');
    end
end
[tmesh0,coeff0,tangent0]=CorrectStartPoint(tmesh0,coeff0,tangent0);
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
s.data.sol=OCMATCONT.formatsolution(tmesh0,coeff0,tangent0);
s.msg='This is the first solution of the BVP continuation';

fprintf('first solution found\n');
fprintf('tangent vector to first solution found\n');
OCMATCONT.testvals=[];
OCMATCONT.uservals=[];
if Singularities
    [tfvals,failed]=EvalTestFunc(ActTest,tmesh0,coeff0,tangent0);
    OCMATCONT.testvals(2,:)=tfvals(ActTest);
    if ~isempty(failed)
        error('Evaluation of test functions failed at starting solution.');
    end
end

sout=s;
if SaveIntermediate
    bvpout.tmesh=tmesh0;
    bvpout.coeff=coeff0;
    bvpout.tangent=tangent0;
    failed=OCMATCONT.saveintermediate(sout,bvpout,1);
end
if HitTargetValue
    % WM: calculate all testfunctions at once
    [hitval failed]=EvalHitFunc(0,tmesh0,coeff0,tangent0);
    if ~isempty(hitval)
        OCMATCONT.hitvzero(2,:)=hitval;
        if ~isempty(failed)
            error('Evaluation of test functions failed at starting solution.');
        end
    else
        HitTargetValue=0;
    end
end

% Open dialog box to interrupt continuation manually
UserStop=stoploop('Stop continuation.');
tmesh1=tmesh0;
coeff1=coeff0;
tangent1=tangent0;
contnum=1;
ind=1;
exitflag=0;
passed=0;
firstnonadmissibleflag=true;
OCMATCONT.contnum=contnum;

while contnum <= MaxContStep && ~exitflag && ~UserStop.Stop()
    corrections=1;
    while 1 && ~UserStop.Stop()
        % predict next solution
        coeffpre=coeff1 + stepwidth * tangent1;

        [tmesh2,coeff2,tangent2,iter]=newtoncorrection(tmesh1,coeffpre,tangent1);

        if ~isempty(coeff2) && ((contnum < dir_check_step) || ((tangent1'*tangent2 > dir_check_angle))) % && tangent2(end)>0
            
            % test solution integrity
            [res,passed]=TestSolutionIntegrity(tmesh2,coeff2);
            nonadmissible=0;
            % test admissibility
            if CheckAdmissibility & contnum>=CheckAdmissibility
                [nonadmissible,infoS,labelS]=EvalAdmissibilityFunc(tmesh2,coeff2,tangent2);
            end
            if ~nonadmissible
                failed=[];
                if Singularities
                    [tfvals,failed]=EvalTestFunc(ActTest,tmesh2,coeff2,tangent2);
                    OCMATCONT.testvals(OCMATCONT.atv,:)=tfvals(ActTest);
                end
                if isempty(failed)||~failed
                    break
                end
            end
            ocmatmsg(' \nNon admissible solution detected, reduce stepwidth.\n')
            if firstnonadmissibleflag
                s.index=contnum;
                s.label='NAS';
                s.data.sol=OCMATCONT.formatsolution(tmesh1,coeff1,tangent1);
                s.data.soln=OCMATCONT.formatsolution(tmesh2,coeff2,tangent2);
                s.data.infoS=infoS;
                s.data.labelS=labelS;
                s.msg='This is the first non-admissible solution.';
                sout=[sout; s];
                firstnonadmissibleflag=false;
            elseif nonadmissible
                s.index=contnum;
                s.label='NAS';
                s.data.sol=OCMATCONT.formatsolution(tmesh1,coeff1,tangent1);
                s.data.soln=OCMATCONT.formatsolution(tmesh2,coeff2,tangent2);
                s.data.infoS=infoS;
                s.data.labelS=labelS;
                s.msg='This is a non-admissible solution.';
                sout=[sout; s];
                % last solution is an admissible solution
            end
        end
        % reset mesh specific entries of global variable OCMATCONT to previous mesh tmesh1
        if WorkSpace
            if OCMATCONT.probleminit(tmesh1,coeff1,tangent1)~=0
                ocmaterror('Initializer failed.');
            end
        end
        if stepwidth > stepwidth_min
            stepwidth=max(stepwidth_min, stepwidth*stepwidth_dec_fac);
            corrections=corrections + 1;
        else      % if not then fail
            ocmatmsg('Current step size too small (point %d)\n',contnum);
            exitflag=1;
            if isempty(coeff2)
                tmesh2=tmesh1;
                coeff2=coeff1;
                tangent2=tangent1;
                MeshAdaptation=0;
                res=0;
            else
                if WorkSpace
                    if OCMATCONT.probleminit(tmesh2,coeff2,tangent2)~=0
                        ocmaterror('Initializer failed.');
                    end
                end
            end
            break
        end
    end
    % show information during continuation on command window
    if PrintContStats && ~UserStop.Stop()
        fprintf(1,'\n Continuation step No.: %i\n',contnum);
        fprintf(1,' stepwidth: %g ([%g - %g])\n',stepwidth,stepwidth_min,stepwidth_max);
        fprintf(1,' Newton Iterations: %g\n',iter);
        fprintf(1,' Mesh size: %g\n',length(tmesh2));
        fprintf(1,' Max. Residuum: %g (%g)\n',max(res),OCMATCONT.OPTIONS.meshadaptreltol);
        PrintContinuationProcess(tmesh2,coeff2,tangent2);
    end
    if PlotCont && ~UserStop.Stop()
        PlotContinuationProcess(tmesh2,coeff2,tangent2);
    end


    % Hit target value
    if HitTargetValue && ~UserStop.Stop()
        [hitval,failed]=EvalHitFunc(0,tmesh2,coeff2,tangent2);
        if isempty(failed) || ~failed
            OCMATCONT.hitvzero(OCMATCONT.ahv,:)=hitval;
            OCMATCONT.ahv=3-OCMATCONT.ahv;
            testchanges=sign(OCMATCONT.hitvzero(1,:))~=sign(OCMATCONT.hitvzero(2,:));
            if any(testchanges)
                testidx=find(testchanges);
                [tmesh2,coeff2,tangent2]=FindSolAtTargetPoint(testidx,tmesh1,coeff1,tangent1,tmesh2,coeff2,tangent2);
                if ~isempty(coeff2)
                    fprintf('\n Target value hit.\n');
                    s.index=contnum;
                    s.label='HTV';
                    s.data.sol=OCMATCONT.formatsolution(tmesh2,coeff2,tangent2);
                    [failed,s]=feval(OCMATCONT.process,-1,tmesh2,coeff2,tangent2,s);
                    s.msg='This is the solution at the target value.';
                    sout=[sout; s];
                    %fprintf(' label=%s\n Continuation parameter=%g\n', s.label,coeff2(OCMATCONT.HE.contparametercoord));
                    if ExitOnTargetValue & contnum>=ExitOnTargetValue
                        exitflag=1;
                    end
                    if PlotCont
                        PlotContinuationProcess(tmesh2,coeff2,tangent2);
                    end

                else
                    tmesh2=tmesh1;
                    coeff2=coeff1;
                    tangent2=tangent1;
                    if ExitOnTargetValue && contnum>=ExitOnTargetValue
                        exitflag=1;
                    end
                end
            end
        end
    end
    
    if stepwidth < stepwidth_max && corrections==1 && iter < OCMATCONT.OPTIONS.maxnewtiter-1
        stepwidth=min(stepwidth*stepwidth_inc_fac, stepwidth_max);
    end
    % Singularities
    %
    if Singularities && ~UserStop.Stop()
        % WM: the testvals arrays are not copied anymore, instead
        % after every iteration the function of both is swapped
        OCMATCONT.atv=3-OCMATCONT.atv;
        % WM: use sign function and compare instead of multiply (for speed).
        testchanges=(sign(OCMATCONT.testvals(1,:)) ~=sign(OCMATCONT.testvals(2,:)));
        if any(testchanges)
            testidx=[ActTest(find( testchanges )) 0]';
            % check if sing occured
            % WM: detect all singularities with a single ismember() call
            %         stz=ismember(cds.SZ,testidx);
            %           singsdetected(ActSing)=all(stz(:,ActSing));
            % DSB: Singularity is detected if
            %  - Every crossing that is required occurs
            %  - Every crossing that is not required does not occur
            %
            % unary plus conversion to double
            S_true=+(OCMATCONT.S(:,ActTest)'==0);  % Required crossings matrix
            S_false=+(OCMATCONT.S(:,ActTest)'==1);  % Required noncrossings matrix
            all_sings_detected=(testchanges*S_true==sum(S_true))&(~testchanges*S_false==sum(S_false));
            singsdetected(ActSing)=all_sings_detected(ActSing);
            if any(singsdetected)
                % singularity detected!
                singsdetected=find(singsdetected==1);
                OCMATCONT.testzero=zeros(OCMATCONT.HE.numdvariables,OCMATCONT.nTest);
                OCMATCONT.testvzero=zeros(OCMATCONT.HE.numdvariables,OCMATCONT.nTest);
                % locate zeros of all testf which changed sign
                %
                tmeshss=[];  % tmesh of singularites
                coeffss=[] ; % v of idem
                tangentss=[] ; % v of idem
                testfound=[];% indices of found zeros of test functions

                sid=[];  % id of idem
                for si=singsdetected
                    %         debug('Singularity %d detected ... ', si);
                    if ismember(si, find(Locators==1))   % do we have a locator?
                        [tmeshs,coeffs,tangents]=OCMATCONT.locate(si,tmesh1,coeff1,tangent1,tmesh2,coeff2,tangent2);
                        lit=0;
                        if ~isempty(coeffs) &&(norm(coeffs-(coeff1+coeff2)/2)<2*norm(coeff2-coeff1))
                            coeffss=[coeffss coeffs];
                            tangentss=[tangentss tangents];
                            sid=[sid si];
                        end
                    else
                        % locate zeros of test functions if not already computed
                        for ti=1:OCMATCONT.nTest
                            if OCMATCONT.S(si,ti)==0 && ~ismember(ti,testfound)
                                [tmeshtf,coefftf,tangenttf,lit]=LocateTestFunction(ti,tmesh1,coeff1,tangent1,tmesh2,coeff2,tangent2);
                                if ~isempty(tmeshtf)
                                    try
                                        d=(norm(coefftf-(coeff1+coeff2)/2)<2*norm(coeff2-coeff1));
                                    catch
                                        [tmesh1tmp coeff1tmp tangent1tmp]=adaptgrid(tmesh1,coeff1,tangent1,tmeshtf);
                                        [tmesh2tmp coeff2tmp tangent2tmp]=adaptgrid(tmesh2,coeff2,tangent2,tmeshtf);
                                        d=(norm(coefftf-(coeff1tmp+coeff2tmp)/2)<2*norm(coeff2tmp-coeff1tmp));
                                    end
                                    if d
                                        %OCMATCONT.testzerosol(ti)=OCMATCONT.formatsolution(tmeshtf,coefftf,tangenttf);
                                        OCMATCONT.testzerosol(ti).tmesh=tmeshtf;
                                        OCMATCONT.testzerosol(ti).coeff=coefftf;
                                        OCMATCONT.testzerosol(ti).tangent=tangenttf;
                                        testfound=[testfound ti];
                                    end
                                    %                               else
                                    %                                   debug('A testfunction for Singularity %d failed to converge ... \n', si);
                                end
                            end
                        end
                        % now we have all zeros/testfunctions we need for
                        % this singularity
                        if any(ismember(testfound, find(OCMATCONT.S(si,:)==0)))
                            [tmeshs,coeffs,tangents]=LocateSingularity(si);
                        else
                            tmeshs=[];
                            coeffs=[];
                            tangents=[];
                        end
                        if ~isempty(tmeshs)
                            tmeshss=[tmeshss;tmeshs];
                            coeffss=[coeffss coeffs];
                            tangentss=[tangentss tangents];
                            sid=[sid si];
                        end
                    end
                end %end of detect/locate loop
                if ~isempty(sid)         % sort
                    [tmeshss,coeffss,tangentss,sid]=xssort(tmesh1,tmeshss,coeffss,tangentss,sid);
                    % WM: moved out of loop for speed
                    sids=1:length(sid);
                    isids=contnum+sids;
                    for si=sids
                        contnum=contnum+1; ind=[ind contnum];s=[];
                        s.msg='';
                        s.index=contnum;
                        s.label=SingLables(sid(si),:);
                        [failed,sf,s]=DefaultProcessor(tmeshss(si,:),coeffss(:,si),tangentss(:,si), s);
                        [tfvals,failed]=EvalTestFunc(ActTest,tmeshss(si,:),coeffss(:,si),tangentss(:,si));
                        s.data.testfunctions=tfvals(ActTest);
                        s.data.sol=OCMATCONT.formatsolution(tmeshss(si,:),coeffss(:,si),tangentss(:,si));
                        [failed,s]=feval(OCMATCONT.process,sid(si),tmeshss(si,:),coeffss(:,si),tangentss(:,si),s);
                        if Userfunctions
                            [ufvals,failed]=feval(OCMATCONT.userfunc, UserInfo, 1:cds.nUserf, tmeshss(si,:), vss(:,si));
                            s.data.userfunctions=ufvals;
                            %                           hout(:,i)=[cds.h;0;s.data.userfunctions';s.data.testfunctions'];
                            hout(:,i)=[0;lit;s.data.userfunctions';s.data.testfunctions'];
                        else
                            % XXX lit=# localisationsteps XXX
                            %                           hout(:,i)=[cds.h;lit;s.data.testfunctions'];
                            %hout(:,i)=[0;lit;s.data.testfunctions'];
                        end
                        sout=[sout; s];
                        if SaveIntermediate
                            bvpout(contnum).tmesh=tmeshss(si,:);
                            bvpout(contnum).coeff=coeffss(:,si);
                            bvpout(contnum).tangent=tangentss(:,si);
                            failed=OCMATCONT.saveintermediate(sout,bvpout,contnum);
                        end
                    end
                end % end of loop over singularities
            end
        end
    end

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    newmesh=false;
    if ((mod(contnum,MeshAdaptation)==0) || ~passed) && ~UserStop.Stop()
        [tmesh3,coeff3,tangent3,newmesh]=AdaptMesh(tmesh2,coeff2,tangent2);
        if length(tmesh3)>OCMATCONT.OPTIONS.maxgridnum
            ocmatmsg([ 'Unable to meet the tolerance without using more than %d '...
                'mesh points. '],OCMATCONT.OPTIONS.maxgridnum);
            exitflag=true;
        else
            failed=DataAdaptation(tmesh3,coeff3,tangent3,tmesh2);
            if ~failed
                tmesh2=tmesh3;
                coeff2=coeff3;
                tangent2=tangent3;
            end
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (mod(contnum,Adapt)==0)
        [res,tmesh2,coeff2,tangent2]=OCMATCONT.adapt(tmesh2,coeff2,tangent2);
        [failed,f]=DefaultProcessor(tmesh2,coeff2,tangent2);
        if res==1 && Singularities
            % recompute testvals
            [tfvals,failed] = EvalTestFunc(ActTest,tmesh2,coeff2,tangent2);
            OCMATCONT.testvals(3-OCMATCONT.atv,:) = tfvals(ActTest);
        end
    end

    if SaveIntermediate
        bvpout(contnum+1).tmesh=tmesh2;
        bvpout(contnum+1).coeff=coeff2;
        bvpout(contnum+1).tangent=tangent2;
        bvpout(contnum+1).stepwidth=stepwidth;
        failed=OCMATCONT.saveintermediate(sout,bvpout,contnum+1);
    end
    contnum=contnum+1;
    OCMATCONT.contnum=contnum;
    % shift tmesh1,tangent1  %
    tmesh1=tmesh2;
    coeff1=coeff2;
    tangent1=tangent2;
end

% The last solution shall pass the integrity test 
if 0%MeshAdaptation
    if newmesh
        while 1
            [tmesh2,coeff2,tangent2]=newtoncorrection(tmesh1,coeff1,tangent1);
            [res,passed]=TestSolutionIntegrity(tmesh2,coeff2);
            if ~passed
                [tmesh1,coeff1,tangent1]=AdaptMesh(tmesh2,coeff2,tangent2);
                failed=DataAdaptation(tmesh1,coeff1,tangent1,tmesh2);
            else
                break
            end
        end
    end
    if length(tmesh1)>OCMATCONT.OPTIONS.maxgridnum
        ocmatmsg([ 'Unable to meet the tolerance without using more than %d '...
            'mesh points. '],OCMATCONT.OPTIONS.maxgridnum);
        exitflag=true;
    else
        tmesh1=tmesh2;
        coeff1=coeff2;
        tangent1=tangent2;
    end
end
fprintf('\n');
EndTime = clock;

DataAdaptation(tmesh1,coeff1,tangent1,tmesh2);
s.index=contnum;
s.label='99';
s.data.sol=OCMATCONT.formatsolution(tmesh1,coeff1,tangent1);
s.msg='This is the last solution of the BVP continuation';
sout=[sout; s];
if SaveIntermediate
    bvpout(contnum+1).tmesh=tmesh2;
    bvpout(contnum+1).coeff=coeff2;
    bvpout(contnum+1).tangent=tangent2;
    failed=OCMATCONT.saveintermediate(sout,bvpout,contnum+1);
end

UserStop.Clear();
clear UserStop
fprintf('elapsed time  = %.1f secs\n', etime(EndTime, StartTime));
%-------------------------------------
%
% Start point Corrector
%
%-------------------------------------

function [tmesh,coeff,tangent]=CorrectStartPoint(tmesh0,coeff0,tangent0)
coeff=[];
tangent=[];
% no tangent vector given, cycle through base-vectors
ndim=length(coeff0);
ii=ndim;
if ~isempty(tangent0)
    [tmesh,coeff,tangent]=newtoncorrection(tmesh0,coeff0,tangent0);
end
tangent0=zeros(ndim,1);
while isempty(coeff) && ii>=1
    tangent0(ii)=1;
    try
        [tmesh,coeff,tangent]=newtoncorrection(tmesh0,coeff0,tangent0);
        if isempty(tmesh)
            return
        end
    catch
        lasterr
    end
    tangent0(ii)=0;
    ii=ii-1;
end

%--< END OF ST PNT CORRECTOR >--


%-------------------------------------
%
% Find the solution at the target point
%
%-------------------------------------
function [tmesh,coeff,tangent]=FindSolAtTargetPoint(id,tmesh1,coeff1,tangent1,tmesh2,coeff2,tangent2)
global OCMATCONT


t1=EvalHitFunc(id,tmesh1,coeff1,tangent1);
t2=EvalHitFunc(id,tmesh2,coeff2,tangent2);

ii=1;
tmesh=[];
tangent=[];
coeff=[];
tmax=10*max(abs(t1(id)),abs(t2(id)));
p=1;
while ii<=OCMATCONT.OPTIONS.maxtestiters
    % WM: make educated guess of where the zero point might be
    if tmax < Inf
        r=abs(t1(id)/(t1(id)-t2(id)))^p;
    else
        r=0.5;
    end
    coeff3=coeff1+r*(coeff2-coeff1);
    tangent3=tangent1+r*(tangent2-tangent1);
    tmesh3=tmesh1;
    [tmesh,coeff,tangent]=newtoncorrection(tmesh3,coeff3,tangent3);

    if isempty(coeff)
        tmesh=tmesh3;
        coeff=coeff3;
        tangent=tangent3;
    end
    tval=EvalHitFunc(id,tmesh,coeff,tangent);
    dist1=norm(coeff-coeff1);
    dist2=norm(coeff-coeff2);
    if abs(tval(id))>tmax
        %     fprintf('testfunction behaving badly.\n');
        tmesh=[];
        tangent=[];
        coeff=[];
        break;
    end
    if abs(tval(id))<=OCMATCONT.OPTIONS.testtolerance && min(dist1,dist2) < 100*OCMATCONT.OPTIONS.vartolerance
        break;
    elseif sign(tval(id))==sign(t2(id))
        coeff2=coeff;
        tangent2=tangent;
        tmesh2=tmesh;
        tmesh1=tmesh;
        t2(id)=tval(id);
        p=1.02;
    else
        coeff1=coeff;
        tangent1=tangent;
        tmesh2=tmesh;
        tmesh1=tmesh;
        t1(id)=tval(id);
        p=0.98;
    end
    ii=ii+1;
    tmesh=[];
    coeff=[];
    tangent=[];
end

%----------------------------------------------
function [tmesh,coeff,tangent,ii]=LocateTestFunction(id,tmesh1,coeff1,tangent1,tmesh2,coeff2,tangent2)
% default locator: bisection
global OCMATCONT


ii=1;
t1=EvalTestFunc(id,tmesh1,coeff1,tangent1);
t2=EvalTestFunc(id,tmesh2,coeff2,tangent2);
tmax=10*max(abs(t1(id)),abs(t2(id)));
p=1;
while ii<=OCMATCONT.OPTIONS.maxtestiters
    % WM: make educated guess of where the zero point might be
    if tmax < Inf
        r=abs(t1(id)/(t1(id)-t2(id)))^p;
    else
        r=0.5;
    end
    coeff3=coeff1+r*(coeff2-coeff1);
    tangent3=tangent1+r*(tangent2-tangent1);
    tmesh3=tmesh1;
    [tmesh,coeff,tangent]=newtoncorrection(tmesh3,coeff3,tangent3);

    if isempty(coeff)
        tmesh=tmesh3;
        coeff=coeff3;
        tangent=tangent3;
    else
        tmesh=tmesh3;
    end
    DefaultProcessor(tmesh,coeff,tangent);
    tval=EvalTestFunc(id,tmesh,coeff,tangent);
    dist1=norm(coeff-coeff1);
    dist2=norm(coeff-coeff2);
    if abs(tval(id)) > tmax
        %     fprintf('testfunction behaving badly.\n');
        tmesh=[];
        tangent=[];
        coeff=[];
        break;
    end
    if abs(tval(id))<=OCMATCONT.OPTIONS.testtolerance && min(dist1,dist2) < 100*OCMATCONT.OPTIONS.vartolerance
        break;
    elseif sign(tval(id))==sign(t2(id))
        coeff2=coeff;
        tangent2=tangent;
        tmesh2=tmesh;
        tmesh1=tmesh;
        t2(id)=tval(id);
        p=1.02;
    else
        coeff1=coeff;
        tangent1=tangent;
        tmesh2=tmesh;
        tmesh1=tmesh;
        t1(id)=tval(id);
        p=0.98;
    end
    ii=ii+1;
    tmesh=[];
    coeff=[];
    tangent=[];
end
%--< END OF locatetestfunction>--
%---------------------------------------------
%----------------------------------------------

%------------------------------------------
%
%  Evaluate testfunctions
%
%------------------------------------------

function [out,failed]=EvalTestFunc(id,tmesh,coeff,tangent)
global OCMATCONT

if id==0
    % WM: evaluate all testfunctions at once
    [out,failed]=feval(OCMATCONT.testfunc,1:OCMATCONT.nTest,tmesh,coeff,tangent);
else
    [out,failed]=feval(OCMATCONT.testfunc,id,tmesh,coeff,tangent);
end

%------------------------------------------
%
%  Evaluate targetvaluefunctions
%
%------------------------------------------

function [out,failed]=EvalHitFunc(id,tmesh,coeff,tangent)
global OCMATCONT

if id==0
    [out,failed]=OCMATCONT.targetvaluefunc(1:OCMATCONT.TargetValueNum,tmesh,coeff,tangent);
else
    [out,failed]=OCMATCONT.targetvaluefunc(id,tmesh,coeff,tangent);
end


%----------------------------
%
% Check admissibility of solution
%
%----------------------------

function [failed,infoS,labelS]=EvalAdmissibilityFunc(tmesh,coeff,tangent)
global OCMATCONT

[failed,infoS,labelS]=OCMATCONT.testadmissibility(tmesh,coeff,tangent);

%--< END OF EVALADMISSIBILITYFUNC >--


%---------------------------------
%
% Plot the solution
%
%
function b=PlotContinuationProcess(tmesh,coeff,tangent)
global OCMATCONT
try
    b=OCMATCONT.plotcontinuation(tmesh,coeff,tangent);
catch
    b=[];
end

%---------------------------------
%
% Print the solution
%
%
function b=PrintContinuationProcess(tmesh,coeff,tangent)
global OCMATCONT
b=OCMATCONT.printcontinuation(tmesh,coeff,tangent);

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

%------------------------------------------------------------
%
% Locate singularity tmeshs between x1 and x2
% First locating zeros of testfunctions, then glue them
%
%------------------------------------------------------------

function [tmeshs,coeffs,tangents]=LocateSingularity(si)
global OCMATCONT

% zero locations of testf i is stored in testzero(:,i)

% if 1 zero/tf check if nonzero/tf _is_ nonzero at that point
% if more zero/tf then glue, nonzero/tf is kind of a problem because can always be found
idx=find( OCMATCONT.S(si,:)==0 );
nzs=find( OCMATCONT.S(si,:)==1 );

len=length(idx);
lnz=length(nzs);

switch len
    case 0
        % Oops, we have detected a singularity without a vanishing testfunction
        error('Internal error: trying to locate a non-detected singularity');

    case 1
        % check all nonzero/tf
        tmeshs=OCMATCONT.testzerosol(idx).tmesh;
        coeffs=OCMATCONT.testzerosol(idx).coeff;
        tangents=OCMATCONT.testzerosol(idx).tangent;

    otherwise
        tmeshz=zeros(len,OCMATCONT.ndim);
        coeffz=zeros(OCMATCONT.ndim,len);
        tangentz=zeros(OCMATCONT.ndim,len);
        nm=zeros(1,len);

        for ii=1:len
            tmeshz(ii,:)=OCMATCONT.testzerosol(idx(ii)).tmesh;
            coeffz(:,ii)=OCMATCONT.testzerosol(idx(ii)).coeff;
            tangentz(:,ii)=OCMATCONT.testzerosol(idx(ii)).tangent;
            nm(ii)=norm(tmeshz(:,ii));
        end

        if max(nm)-min(nm) < OCMATCONT.OPTIONS.funtolerance
            tmeshs=mean(tmeshz,1)';
            coeffs=mean(coeffz',1)';
            tangents=mean(tangentz',1)';
        else
            tmeshs=[];
            coeffs=[];
            tangents=[];
            return;
        end
end

if lnz==0, return; end

% checking non zeros

DefaultProcessor(tmeshs,coeffs,tangents);
tval=EvalTestFunc(nzs,tmeshs,coeffs,tangents);
if any(abs(tval(nzs))<=OCMATCONT.OPTIONS.testtolerance)
    tmeshs=[];
    coeffs=[];
    tangents=[];
end


%------------------------------------------------
%
%  Sorts [xs,vs] and sing with x1 starting point
%  xs contains x_singular, vs v_singular
%  sing contains their id's
%  Sort criterion: dist(x1,x)
%
%------------------------------------------------

function [tmeshs,coeffs,tangents,sing]=xssort(tmesh1,tmeshs,coeffs,tangents,sing)
% WM: Matlab has a sort function, beter use it...
len=size(tmeshs,1);
if len > 1
    tmesho=tmesh1(ones(1,len),:);
    [dummy,i]=sortrows((tmeshs-tmesho));
    tmeshs=tmeshs(i,:);
    coeffs=coeffs(:,i);
    tangents=tangents(:,i);
    sing=sing(:,i);
end

%--< END OF XSSORT >--

%----------------------------
%
% DefaultProcessor
%
%----------------------------

function [failed,f,s]=DefaultProcessor(tmesh,coeff,tangent,s)
global OCMATCONT
% WM: this now actually calls the default processor,
% either with or without a singular point structure

if nargin > 3
    [failed,f,s]=OCMATCONT.defaultprocessor(tmesh,coeff,tangent,s);
else
    [failed,f]=OCMATCONT.defaultprocessor(tmesh,coeff,tangent);
end


%------------------------------------------------
%
% Call method specific mesh adaptation
%
%------------------------------------------------


function [tmesh,coeff,tangent,newmesh]=AdaptMesh(tmesh,coeff,tangent)
global OCMATCONT
newmesh=false;
[t,y,z,freepar,modelpar]=OCMATCONT.drearr(tmesh,coeff);
switch OCMATCONT.bvpmethod
    case 'bvp6c'
        res=residual_bvp6c(t,y,freepar,modelpar,OCMATCONT.RHS(1:OCMATCONT.HE.numdvariablesmc),OCMATCONT.ode);
        if max(res) < OCMATCONT.OPTIONS.meshadaptreltol;
            return
        end
    case 'bvp5c'
        res=residual_bvp5c(t,y,freepar,modelpar,OCMATCONT.RHS(1:OCMATCONT.HE.numdvariablesmc),OCMATCONT.ode);
        if verifyresidual_bvp5c(max(res))
            return
        end
    case 'bvp4c'
        if isempty(OCMATCONT.RHS)
            OCMATCONT.RHS=OCMATCONT.operatoreq(tmesh,coeff,[],OCMATCONT.ode,OCMATCONT.bc,OCMATCONT.icfun);
        end
        res=residual_bvp4c(t,y,freepar,modelpar,OCMATCONT.RHS(1:OCMATCONT.HE.numdvariablesmc),OCMATCONT.ode);
        if max(res) < OCMATCONT.OPTIONS.meshadaptreltol;
            [tmesh,ynew,tangent]=removepoints_bvp4c(t,y,tangent,res);
            coeff=[ynew(:);freepar(:)];
            return
        end
    case 'gbvp4c'
        if isempty(OCMATCONT.RHS)
            OCMATCONT.RHS=OCMATCONT.operatoreq(tmesh,coeff,[],OCMATCONT.ode,OCMATCONT.bc,OCMATCONT.icfun);
        end
        res=residual_gbvp4c(t,y,freepar,modelpar,OCMATCONT.RHS(1:OCMATCONT.HE.numdvariablesmc),OCMATCONT.ode);
        if max(res) < OCMATCONT.OPTIONS.meshadaptreltol;
            [tmesh,ynew,tangent]=removepoints_gbvp4c(t,y,tangent,res);
            dataadaptation_gbvp4c(tmesh)
            coeff=ynew(OCMATCONT.HE.DDATA.meshvalcoord);
            coeff=[coeff(:);freepar(:)];
            return
        end
end
% Detect mesh oscillations:  Was there a mesh with
% the same number of nodes and a similar residual?
% residualReduction = abs(OCMATCONT.meshHistory((OCMATCONT.meshHistory(:,1) == OCMATCONT.HE.TIMEDDATA.nummesh),2) - maxres)/maxres;
oscLikely = false;%any( residualReduction < OCMATCONT.OPTIONS.residualreductionguard);

switch OCMATCONT.bvpmethod
    case 'bvp6c'
        [tmesh,ynew,tangent]=meshadaptation_bvp6c(t,y,tangent,res,~oscLikely);
        coeff=[ynew(:);freepar(:)];
    case 'bvp5c'
        [tmesh,ynew,tangent]=meshadaptation_bvp5c(t,y,tangent,res,~oscLikely,freepar,modelpar,OCMATCONT.ode);
        coeff=[ynew(:)];
    case 'bvp4c'
        [tmesh,ynew,tangent]=meshadaptation_bvp4c(t,y,tangent,res,~oscLikely);
        coeff=[ynew(:);freepar(:)];
    case 'gbvp4c'
        [tmesh,ynew,tangent]=meshadaptation_gbvp4c(t,y,tangent,res,~oscLikely);
        dataadaptation_gbvp4c(tmesh)
        coeff=ynew(OCMATCONT.HE.DDATA.meshvalcoord);
        coeff=[coeff(:);freepar(:)];
end
% if ~isempty(OCMATCONT.dataadaptation)
%     OCMATCONT.dataadaptation(tmesh,t);
% end
newmesh=true;



function [res,passed]=TestSolutionIntegrity(tmesh,coeff)
global OCMATCONT
[t,y,z,freepar,modelpar]=OCMATCONT.drearr(tmesh,coeff);
switch OCMATCONT.bvpmethod
    case 'bvp6c'
        res=residual_bvp6c(t,y,freepar,modelpar,OCMATCONT.RHS(1:OCMATCONT.HE.numdvariablesmc),OCMATCONT.ode);
        passed=max(res) < OCMATCONT.OPTIONS.meshadaptreltol;
    case 'bvp5c'
        res=residual_bvp5c(t,y,freepar,modelpar,OCMATCONT.RHS(1:OCMATCONT.HE.numdvariablesmc),OCMATCONT.ode);
        passed=verifyresidual_bvp5c(max(res));
    case 'bvp4c'
        res=residual_bvp4c(t,y,freepar,modelpar,OCMATCONT.RHS(1:OCMATCONT.HE.numdvariablesmc),OCMATCONT.ode);
        passed=max(res) < OCMATCONT.OPTIONS.meshadaptreltol;
    case 'gbvp4c'
        res=residual_gbvp4c(t,y,freepar,modelpar,OCMATCONT.RHS(1:OCMATCONT.HE.numdvariablesmc),OCMATCONT.ode);
        passed=max(res) < OCMATCONT.OPTIONS.meshadaptreltol;
end

