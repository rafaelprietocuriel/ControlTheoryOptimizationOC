function sout=daecont(varargin)
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
OCMATCONT.settargetvalue=problemhandles{10};
OCMATCONT.operatorpfrechet=problemhandles{11};
OCMATCONT.residual=problemhandles{12};
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
OCMATCONT.writelogfile=problemhandles{27};

dae=problemhandles{4}{1};
bc=problemhandles{4}{2};
ic=problemhandles{4}{3};
if getocoptions(opt,'SBVPOC','FJacobian')
    daejac=problemhandles{5}{1};
else
    daejac=[];
end
if getocoptions(opt,'SBVPOC','BCJacobian')
    bcjac=problemhandles{5}{2};
else
    bcjac=[];
end
daehess=[];
bchess=[];
daetensor3=[];
icfunjac=[];
bctensor3=[];

OCMATCONT.OPTIONS.newtonabstol=getocoptions(opt,'NEWTON','AbsTol');
OCMATCONT.OPTIONS.newtonreltol=getocoptions(opt,'NEWTON','RelTol');
OCMATCONT.OPTIONS.maxnewtiter=getocoptions(opt,'NEWTON','MaxNewtonIters');
OCMATCONT.OPTIONS.switchtoffnfactor=getocoptions(opt,'NEWTON','SwitchToFFNFactor');
OCMATCONT.OPTIONS.lambdamin=getocoptions(opt,'NEWTON','LambdaMin');
OCMATCONT.OPTIONS.maxprobes=getocoptions(opt,'NEWTON','MaxProbes');
OCMATCONT.OPTIONS.updatejacfactor=getocoptions(opt,'NEWTON','UpdateJacFactor');
OCMATCONT.OPTIONS.checksingular=logical(getocoptions(opt,'NEWTON','CheckSingular'));
OCMATCONT.OPTIONS.display=getocoptions(opt,'NEWTON','Display');
%OCMATCONT.OPTIONS.singularthreshold=getocoptions(opt,'NEWTON','SingularThreshold');
OCMATCONT.newtonsolver=str2func(getocoptions(opt,'GENERAL','NewtonSolver'));
OCMATCONT.OPTIONS.admissibletol=getocoptions(opt,'GENERAL','AdmissibleTolerance');
OCMATCONT.OPTIONS.zerotimedifftolerance=getocoptions(opt,'GENERAL','ZeroTimeDiffTolerance');
OCMATCONT.OPTIONS.contlog=getocoptions(opt,'GENERAL','ContLog');

OCMATCONT.OPTIONS.meshadaptabstol=getocoptions(opt,'SBVPOC','MeshAdaptAbsTol');
OCMATCONT.OPTIONS.meshadaptreltol=getocoptions(opt,'SBVPOC','MeshAdaptRelTol');
OCMATCONT.OPTIONS.meshadaptfinemesh=getocoptions(opt,'SBVPOC','MeshAdaptFineMesh');
OCMATCONT.OPTIONS.meshupdatemode=getocoptions(opt,'SBVPOC','MeshUpdateMode');
OCMATCONT.OPTIONS.meshadaptfactormax=getocoptions(opt,'SBVPOC','MeshAdaptFactorMax');
OCMATCONT.OPTIONS.meshadaptmaxiter=getocoptions(opt,'SBVPOC','MeshAdaptMaxIter');
OCMATCONT.OPTIONS.maxgridnum=getocoptions(opt,'SBVPOC','NMax');
OCMATCONT.OPTIONS.xyvectorized=strcmpi(getocoptions(opt,'BVP','Vectorized'),'on');
OCMATCONT.OPTIONS.krange=getocoptions(opt,'SBVPOC','KRange');

OCMATCONT.OPTIONS.totalrelativedistance=getocoptions(opt,'OCCONTARG','TotalRelativeDistance');
OCMATCONT.OPTIONS.maxtestiters=getocoptions(opt,'OCCONTARG','MaxTestIters');
OCMATCONT.OPTIONS.testtolerance=getocoptions(opt,'OCCONTARG','TestTolerance');
OCMATCONT.OPTIONS.funtolerance=getocoptions(opt,'OCCONTARG','FunTolerance');
OCMATCONT.OPTIONS.increment=getocoptions(opt,'OCCONTARG','Increment');
OCMATCONT.OPTIONS.vartolerance=getocoptions(opt,'OCCONTARG','VarTolerance');
OCMATCONT.OPTIONS.saveintermediate=strcmp(getocoptions(opt,'OCCONTARG','SaveIntermediate'),'on');
OCMATCONT.OPTIONS.continuationmethod=getocoptions(opt,'OCCONTARG','ContinuationMethod'); % 0 ... Arc-length continuation, 1 Moore-Penrose continuation
OCMATCONT.OPTIONS.maxdistance=getocoptions(opt,'OCCONTARG','MaxDistance');

OCMATCONT.bvpmethod='sbvpoc';%getocoptions(opt,'GENERAL','BVPMethod');

% Handle argument functions and additional arguments

OCMATCONT.dae=dae;
OCMATCONT.bc=bc;
OCMATCONT.icfun=ic;
OCMATCONT.daejac=daejac;
OCMATCONT.bcjac=bcjac;
OCMATCONT.icfunjac=icfunjac;
OCMATCONT.daehess=daehess;
OCMATCONT.bchess=bchess;
OCMATCONT.daetensor3=daetensor3;
OCMATCONT.bctensor3=bctensor3;

tcolmesh0=sol0.data.xcol;
coeff0=sol0.data.coeff;

if OCMATCONT.OPTIONS.contlog
    failed=OCMATCONT.writelogfile([]);
    if failed
        OCMATCONT.OPTIONS.contlog=0;
    else
        failed=OCMATCONT.writelogfile(sprintf('Newtonsolver: %s\n',func2str(OCMATCONT.newtonsolver)));
        fn=fieldnames(OCMATCONT.OPTIONS);
        for ii=1:length(fn)
            val=OCMATCONT.OPTIONS.(fn{ii});
            val=val(:).';
            if ~isempty(val) && rem(val(1),1)==0
                OCMATCONT.writelogfile(sprintf('%s:\t',fn{ii}));
                OCMATCONT.writelogfile(sprintf('%d',val));
                OCMATCONT.writelogfile(sprintf('\n'));
            elseif ~isempty(val) && abs(val(1))<1e-4
                OCMATCONT.writelogfile(sprintf('%s:\t',fn{ii}));
                OCMATCONT.writelogfile(sprintf('%4',val));
                OCMATCONT.writelogfile(sprintf('\n'));
            else
                OCMATCONT.writelogfile(sprintf('%s:\t',fn{ii}));
                OCMATCONT.writelogfile(sprintf('%f',val));
                OCMATCONT.writelogfile(sprintf('\n'));
            end
        end
    end
end
if isfield(sol0.data,'tangent')
    tangent0=sol0.data.tangent0;
else
    tangent0=[];
end

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
        OCMATCONT.testzerosol(OCMATCONT.nActTest).tcolmesh=[];
        OCMATCONT.testzerosol(OCMATCONT.nActTest).coeff=[];
        OCMATCONT.testzerosol(OCMATCONT.nActTest).tangent=[];
        if isempty(ActTest)
            Singularities=0;
        end
    end
end
%HitTargetValue=0; % testing
HitTargetValue=HitTargetValue&&OCMATCONT.TargetValueNum>0;
if HitTargetValue
    OCMATCONT.TargetValueNum=1;
    OCMATCONT.hitvzero=zeros(2,OCMATCONT.TargetValueNum);  % tangent where testf is zero
    OCMATCONT.ahv=1;
end

if PlotCont
    clf
end
StartTime = clock;

% initialize user workspace
% if WorkSpace
%     if OCMATCONT.probleminit(tcolmesh0,coeff0,tangent0)~=0
%         ocmaterror('Initializer failed.');
%     end
% end
if OCMATCONT.OPTIONS.contlog
    OCMATCONT.writelogfile(sprintf('Correct Start Point\n'));
end

[tcolmesh0,coeff0,tangent0]=CorrectStartPoint(tcolmesh0,coeff0,tangent0);
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
s.data.sol=OCMATCONT.formatsolution(tcolmesh0,coeff0,tangent0);
s.msg='This is the first solution of the BVP continuation';

fprintf('first solution found\n');
fprintf('tangent vector to first solution found\n');
OCMATCONT.testvals=[];
OCMATCONT.uservals=[];
if Singularities
    [tfvals,failed]=EvalTestFunc(ActTest,tcolmesh0,coeff0,tangent0);
    OCMATCONT.testvals(2,:)=tfvals(ActTest);
    if ~isempty(failed)
        error('Evaluation of test functions failed at starting solution.');
    end
end

sout=s;
if SaveIntermediate
    bvpout.tcolmesh=tcolmesh0;
    bvpout.coeff=coeff0;
    bvpout.tangent=tangent0;
    failed=OCMATCONT.saveintermediate(sout,bvpout,1);
end
if HitTargetValue
    % WM: calculate all testfunctions at once
    [hitval,failed]=EvalHitFunc(0,tcolmesh0,coeff0,tangent0);
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
tcolmesh1=tcolmesh0;
coeff1=coeff0;
tangent1=tangent0;
contnum=1;
ind=1;
exitflag=0;
OCMATCONT.contnum=contnum;

while contnum <= MaxContStep && ~exitflag && ~UserStop.Stop()
    corrections=1;
    while 1 && ~UserStop.Stop()
        if OCMATCONT.OPTIONS.contlog
            failed=OCMATCONT.writelogfile(sprintf('Continuation Step: %d,\t\tStepwidth: %f\n',contnum,stepwidth));
        end
        % predict next solution
        coeffpre=coeff1+stepwidth*tangent1;

        [tcolmesh2,coeff2,tangent2,iter]=newtoncorrection(tcolmesh1,coeffpre,tangent1);

        if ~isempty(coeff2) && ((contnum < dir_check_step) || ((tangent1'*tangent2 > dir_check_angle))) % && tangent2(end)>0

            % test solution integrity
            if ~isnan(mod(contnum,MeshAdaptation)) && ~mod(contnum,MeshAdaptation)
                if PrintContStats
                    fprintf(1,'\n Start mesh adaptation:\n');
                end
                [tcolmesh3,coeff3,tangent3,passed,errnrm,errnrm2]=meshadaptation4dae(tcolmesh2,coeff2,tangent2);
                if ~isempty(tcolmesh3)
                    tcolmesh2=tcolmesh3;
                    coeff2=coeff3;
                    tangent2=tangent3;
                    if PrintContStats
                        fprintf(1,' Mesh adaptation finished succesfully.\n');
                    end
                    % if passed = 0, reduce stepsize starting with new
                    % solution.
                else
                    tmesh1=tcolmesh1(1:OCMATCONT.CollocationNumber+1:end);
                    updatedaecoefficients(tmesh1,'MESHDATA');
                    tmesh11_2=makefinemesh(tmesh1);
                    updatedaecoefficients(tmesh11_2,'MESHDATA1_2');
                    if PrintContStats
                        fprintf(1,'Mesh adaptation was not succesfull.\n');
                    end
                end
            else
                errnrm=[];
                passed=1;
            end
            if passed
                break
            end
        end
        % reset mesh specific entries of global variable OCMATCONT to previous mesh tcolmesh1
        %         if WorkSpace
        %             if OCMATCONT.probleminit(tcolmesh1,coeff1,tangent1)~=0
        %                 ocmaterror('Initializer failed.');
        %             end
        %         end
        if PrintContStats
            fprintf(1,'Step size reduced.\n');
        end
        if stepwidth > stepwidth_min
            stepwidth=max(stepwidth_min, stepwidth*stepwidth_dec_fac);
            corrections=corrections + 1;
        else      % if not then fail
            ocmatmsg('Current step size too small (point %d)\n',contnum);
            exitflag=1;
            if isempty(coeff2)
                tcolmesh2=tcolmesh1;
                coeff2=coeff1;
                tangent2=tangent1;
                MeshAdaptation=0;
            else
%                 if WorkSpace
%                     if OCMATCONT.probleminit(tcolmesh2,coeff2,tangent2)~=0
%                         ocmaterror('Initializer failed.');
%                     end
%                 end
            end
            break
        end
    end
    % show information during continuation on command window
    if PrintContStats
        fprintf(1,'\n Continuation step No.: %i\n',contnum);
        fprintf(1,' stepwidth: %g ([%g - %g])\n',stepwidth,stepwidth_min,stepwidth_max);
        fprintf(1,' Newton Iterations: %g (%g)\n',iter,OCMATCONT.OPTIONS.maxnewtiter);
        fprintf(1,' Mesh size: %g\n',OCMATCONT.MESHDATA.meshNumber);
        if ~isempty(errnrm)
            if length(OCMATCONT.OPTIONS.meshadaptreltol)>1
                fprintf(1,' Max. Error norm: %g rel(%g-%g) abs(%g-%g)\n',errnrm,min(OCMATCONT.OPTIONS.meshadaptreltol),max(OCMATCONT.OPTIONS.meshadaptreltol),min(OCMATCONT.OPTIONS.meshadaptabstol),max(OCMATCONT.OPTIONS.meshadaptabstol));
                fprintf(1,' Max. Error norm: %g rel(%g-%g) abs(%g-%g)\n',errnrm2,min(OCMATCONT.OPTIONS.meshadaptreltol),max(OCMATCONT.OPTIONS.meshadaptreltol),min(OCMATCONT.OPTIONS.meshadaptabstol),max(OCMATCONT.OPTIONS.meshadaptabstol));
            else
                fprintf(1,' Max. Error norm: %g rel(%g) abs(%g)\n',errnrm,OCMATCONT.OPTIONS.meshadaptreltol,OCMATCONT.OPTIONS.meshadaptabstol);
                fprintf(1,' Max. Error norm: %g rel(%g) abs(%g)\n',errnrm2,OCMATCONT.OPTIONS.meshadaptreltol,OCMATCONT.OPTIONS.meshadaptabstol);
            end
        end
        PrintContinuationProcess(tcolmesh2,coeff2,tangent2);
    end
    if PlotCont
        PlotContinuationProcess(tcolmesh2,coeff2,tangent2);
    end
    % Hit target value
    if HitTargetValue
        [hitval,failed]=EvalHitFunc(0,tcolmesh2,coeff2,tangent2);
        if isempty(failed) || ~failed
            OCMATCONT.hitvzero(OCMATCONT.ahv,:)=hitval;
            OCMATCONT.ahv=3-OCMATCONT.ahv;
            testchanges=sign(OCMATCONT.hitvzero(1,:))~=sign(OCMATCONT.hitvzero(2,:));
            if any(testchanges)
                testidx=find(testchanges);
                [tcolmesh2,coeff2,tangent2]=ApproximateSolAtTargetPoint(testidx,tcolmesh1,coeff1,tangent1,tcolmesh2,coeff2,tangent2);
                if ~isempty(coeff2)
                    [tcolmesh3,coeff3,tangent3]=ExactSolAtTargetPoint(testidx,tcolmesh2,coeff2,tangent2);
                    if ~isempty(tcolmesh3)
                        tcolmesh2=tcolmesh3;
                        coeff2=coeff3;
                        tangent2=tangent3;
                    end
                end
                if ~isempty(coeff2)
                    
                    fprintf('\n Target value hit.\n');
                    s.index=contnum;
                    s.label='HTV';
                    s.data.sol=OCMATCONT.formatsolution(tcolmesh2,coeff2,tangent2);
                    [failed,s]=feval(OCMATCONT.process,-1,tcolmesh2,coeff2,tangent2,s);
                    s.msg='This is the solution at the target value.';
                    sout=[sout; s];
                    %fprintf(' label=%s\n Continuation parameter=%g\n', s.label,coeff2(OCMATCONT.HE.contparametercoord));
                    if ExitOnTargetValue && contnum>=ExitOnTargetValue
                        exitflag=1;
                    end
                    if PlotCont
                        PlotContinuationProcess(tcolmesh2,coeff2,tangent2);
                    end

                else
                    tcolmesh2=tcolmesh1;
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
    if Singularities
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
                tcolmeshss=[];  % tcolmesh of singularites
                coeffss=[] ; % v of idem
                tangentss=[] ; % v of idem
                testfound=[];% indices of found zeros of test functions

                sid=[];  % id of idem
                for si=singsdetected
                    %         debug('Singularity %d detected ... ', si);
                    if ismember(si, find(Locators==1))   % do we have a locator?
                        [tcolmeshs,coeffs,tangents]=OCMATCONT.locate(si,tcolmesh1,coeff1,tangent1,tcolmesh2,coeff2,tangent2);
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
                                [tcolmeshtf,coefftf,tangenttf,lit]=LocateTestFunction(ti,tcolmesh1,coeff1,tangent1,tcolmesh2,coeff2,tangent2);
                                if ~isempty(tcolmeshtf)
                                    try
                                        d=(norm(coefftf-(coeff1+coeff2)/2)<2*norm(coeff2-coeff1));
                                    catch
                                        [tcolmesh1tmp coeff1tmp tangent1tmp]=adaptgrid(tcolmesh1,coeff1,tangent1,tcolmeshtf);
                                        [tcolmesh2tmp coeff2tmp tangent2tmp]=adaptgrid(tcolmesh2,coeff2,tangent2,tcolmeshtf);
                                        d=(norm(coefftf-(coeff1tmp+coeff2tmp)/2)<2*norm(coeff2tmp-coeff1tmp));
                                    end
                                    if d
                                        %OCMATCONT.testzerosol(ti)=OCMATCONT.formatsolution(tcolmeshtf,coefftf,tangenttf);
                                        OCMATCONT.testzerosol(ti).tcolmesh=tcolmeshtf;
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
                            [tcolmeshs,coeffs,tangents]=LocateSingularity(si);
                        else
                            tcolmeshs=[];
                            coeffs=[];
                            tangents=[];
                        end
                        if ~isempty(tcolmeshs)
                            tcolmeshss=[tcolmeshss;tcolmeshs];
                            coeffss=[coeffss coeffs];
                            tangentss=[tangentss tangents];
                            sid=[sid si];
                        end
                    end
                end %end of detect/locate loop
                if ~isempty(sid)         % sort
                    [tcolmeshss,coeffss,tangentss,sid]=xssort(tcolmesh1,tcolmeshss,coeffss,tangentss,sid);
                    % WM: moved out of loop for speed
                    sids=1:length(sid);
                    isids=contnum+sids;
                    for si=sids
                        contnum=contnum+1; ind=[ind contnum];s=[];
                        s.index=contnum;
                        s.label=SingLables(sid(si),:);
                        [failed,sf,s]=DefaultProcessor(tcolmeshss(si,:),coeffss(:,si),tangentss(:,si), s);
                        [tfvals,failed]=EvalTestFunc(ActTest,tcolmeshss(si,:),coeffss(:,si),tangentss(:,si));
                        s.data.testfunctions=tfvals(ActTest);
                        s.data.sol=OCMATCONT.formatsolution(tcolmeshss(si,:),coeffss(:,si),tangentss(:,si));
                        [failed,s]=feval(OCMATCONT.process,sid(si),tcolmeshss(si,:),coeffss(:,si),tangentss(:,si),s);
                        if Userfunctions
                            [ufvals,failed]=feval(OCMATCONT.userfunc, UserInfo, 1:cds.nUserf, tcolmeshss(si,:), vss(:,si));
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
                            bvpout(contnum).tcolmesh=tcolmeshss(si,:);
                            bvpout(contnum).coeff=coeffss(:,si);
                            bvpout(contnum).tangent=tangentss(:,si);
                            failed=OCMATCONT.saveintermediate(sout,bvpout,contnum);
                        end
                    end
                end % end of loop over singularities
            end
        end
    end

    if SaveIntermediate
        bvpout(contnum+1).tcolmesh=tcolmesh2;
        bvpout(contnum+1).coeff=coeff2;
        bvpout(contnum+1).tangent=tangent2;
        bvpout(contnum+1).stepwidth=stepwidth;
        failed=OCMATCONT.saveintermediate(sout,bvpout,contnum+1);
    end
    contnum=contnum+1;
    OCMATCONT.contnum=contnum;
    % shift tcolmesh1,tangent1  %
    tcolmesh1=tcolmesh2;
    coeff1=coeff2;
    tangent1=tangent2;
end
fprintf('\n');
EndTime = clock;

s.index=contnum;
s.label='99';
s.data.sol=OCMATCONT.formatsolution(tcolmesh1,coeff1,tangent1);
s.msg='This is the last solution of the BVP continuation';
sout=[sout; s];
if SaveIntermediate
    bvpout(contnum+1).tcolmesh=tcolmesh2;
    bvpout(contnum+1).coeff=coeff2;
    bvpout(contnum+1).tangent=tangent2;
    failed=OCMATCONT.saveintermediate(sout,bvpout,contnum+1);
end
if OCMATCONT.OPTIONS.contlog
      failed=OCMATCONT.writelogfile(-1);
end
UserStop.Clear();
clear UserStop
fprintf('elapsed time  = %.1f secs\n', etime(EndTime, StartTime));
%-------------------------------------
%
% Start point Corrector
%
%-------------------------------------

function [tcolmesh,coeff,tangent,method]=CorrectStartPoint(tcolmesh0,coeff0,tangent0)
coeff=[];
tangent=[];
% no tangent vector given, cycle through base-vectors
ndim=length(coeff0);
ii=ndim;
if ~isempty(tangent0)
    [tcolmesh,coeff,tangent]=newtoncorrection(tcolmesh0,coeff0,tangent0);
end
tangent0=zeros(ndim,1);
while isempty(coeff) && ii>=1
    tangent0(ii)=1;
    try
        [tcolmesh,coeff,tangent]=newtoncorrection(tcolmesh0,coeff0,tangent0);
    catch
        lasterr
    end
    tangent0(ii)=0;
    ii=ii-1;
end

%--< END OF ST PNT CORRECTOR >--


%-------------------------------------
%
% Approximate the solution at the target point
%
%-------------------------------------
function [tcolmesh,coeff,tangent]=ApproximateSolAtTargetPoint(id,tcolmesh1,coeff1,tangent1,tcolmesh2,coeff2,tangent2)
global OCMATCONT

if length(tcolmesh2)~=length(tcolmesh1)
    MESHDATA=OCMATCONT.MESHDATA;
    MESHDATA1_2=OCMATCONT.MESHDATA1_2;
    tmesh1=tcolmesh1(1:OCMATCONT.CollocationNumber+1:end);
    tmesh2=tcolmesh2(1:OCMATCONT.CollocationNumber+1:end);
    updatedaecoefficients(tmesh1,'MESHDATA');
    tmesh11_2=makefinemesh(tmesh1);
    updatedaecoefficients(tmesh11_2,'MESHDATA1_2');
    sol=sol2daestruct(0,tcolmesh1,coeff1);
    coeff1=points2coeff(sol,tmesh2);
    coeff1(end-OCMATCONT.freeparameternum+1:end)=sol.parameters;

    tcolmesh1=tcolmesh2;
    OCMATCONT.MESHDATA=MESHDATA;
    OCMATCONT.MESHDATA1_2=MESHDATA1_2;
    tangent1=calculatenewtangent(tcolmesh1,coeff1,sign(tangent1(end)));
    [tcolmesh1,coeff1]=newtoncorrection(tcolmesh1,coeff1,[]);
    if isempty(tcolmesh1)
        tcolmesh=[];
        coeff=[];
        tangent=[];
        return
    end
end
t1=EvalHitFunc(id,tcolmesh1,coeff1,tangent1);
t2=EvalHitFunc(id,tcolmesh2,coeff2,tangent2);

ii=1;
tcolmesh=[];
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
    tcolmesh3=tcolmesh1;
    [tcolmesh,coeff,tangent]=newtoncorrection(tcolmesh3,coeff3,tangent3);

    if isempty(coeff)
        tcolmesh=tcolmesh3;
        coeff=coeff3;
        tangent=tangent3;
    end
    tval=EvalHitFunc(id,tcolmesh,coeff,tangent);
    dist1=norm(coeff-coeff1);
    dist2=norm(coeff-coeff2);
    if abs(tval(id))>tmax
        %     fprintf('testfunction behaving badly.\n');
        tcolmesh=[];
        tangent=[];
        coeff=[];
        break;
    end
    if abs(tval(id))<=OCMATCONT.OPTIONS.testtolerance && min(dist1,dist2) < 100*OCMATCONT.OPTIONS.vartolerance
        break;
    elseif sign(tval(id))==sign(t2(id))
        coeff2=coeff;
        tangent2=tangent;
        %tcolmesh2=tcolmesh;
        tcolmesh1=tcolmesh;
        t2(id)=tval(id);
        p=1.02;
    else
        coeff1=coeff;
        tangent1=tangent;
        %tcolmesh2=tcolmesh;
        tcolmesh1=tcolmesh;
        t1(id)=tval(id);
        p=0.98;
    end
    ii=ii+1;
    tcolmesh=[];
    coeff=[];
    tangent=[];
end



%-------------------------------------
%
% Calculate the solution at the target point
%
%-------------------------------------
function [tcolmesh,coeff,tangent]=ExactSolAtTargetPoint(id,tcolmesh,coeff,tangent)
global OCMATCONT

coeff=OCMATCONT.settargetvalue(id,tcolmesh,coeff);
[tcolmesh1,coeff1]=newtoncorrection(tcolmesh,coeff,[]);
if isempty(tcolmesh1)
    tcolmesh=[];
    coeff=[];
    tangent=[];
    return
end
tcolmesh=tcolmesh1;
coeff=coeff1;
tangent=calculatenewtangent(tcolmesh1,coeff1,sign(tangent(end)));

%----------------------------------------------
function [tcolmesh,coeff,tangent,ii]=LocateTestFunction(id,tcolmesh1,coeff1,tangent1,tcolmesh2,coeff2,tangent2)
% default locator: bisection
global OCMATCONT


ii=1;
t1=EvalTestFunc(id,tcolmesh1,coeff1,tangent1);
t2=EvalTestFunc(id,tcolmesh2,coeff2,tangent2);
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
    tcolmesh3=tcolmesh1;
    [tcolmesh,coeff,tangent]=newtoncorrection(tcolmesh3,coeff3,tangent3);

    if isempty(coeff)
        tcolmesh=tcolmesh3;
        coeff=coeff3;
        tangent=tangent3;
    else
        tcolmesh=tcolmesh3;
    end
    DefaultProcessor(tcolmesh,coeff,tangent);
    tval=EvalTestFunc(id,tcolmesh,coeff,tangent);
    dist1=norm(coeff-coeff1);
    dist2=norm(coeff-coeff2);
    if abs(tval(id)) > tmax
        %     fprintf('testfunction behaving badly.\n');
        tcolmesh=[];
        tangent=[];
        coeff=[];
        break;
    end
    if abs(tval(id))<=OCMATCONT.OPTIONS.testtolerance && min(dist1,dist2) < 100*OCMATCONT.OPTIONS.vartolerance
        break;
    elseif sign(tval(id))==sign(t2(id))
        coeff2=coeff;
        tangent2=tangent;
        tcolmesh2=tcolmesh;
        tcolmesh1=tcolmesh;
        t2(id)=tval(id);
        p=1.02;
    else
        coeff1=coeff;
        tangent1=tangent;
        tcolmesh2=tcolmesh;
        tcolmesh1=tcolmesh;
        t1(id)=tval(id);
        p=0.98;
    end
    ii=ii+1;
    tcolmesh=[];
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

function [out,failed]=EvalTestFunc(id,tcolmesh,coeff,tangent)
global OCMATCONT

if id==0
    % WM: evaluate all testfunctions at once
    [out,failed]=feval(OCMATCONT.testfunc,1:OCMATCONT.nTest,tcolmesh,coeff,tangent);
else
    [out,failed]=feval(OCMATCONT.testfunc,id,tcolmesh,coeff,tangent);
end

%------------------------------------------
%
%  Evaluate targetvaluefunctions
%
%------------------------------------------

function [out,failed]=EvalHitFunc(id,tcolmesh,coeff,tangent)
global OCMATCONT

if id==0
    [out,failed]=OCMATCONT.targetvaluefunc(1:OCMATCONT.TargetValueNum,tcolmesh,coeff,tangent);
else
    [out,failed]=OCMATCONT.targetvaluefunc(id,tcolmesh,coeff,tangent);
end


%----------------------------
%
% Check admissibility of solution
%
%----------------------------

function [failed,infoS,labelS]=EvalAdmissibilityFunc(tcolmesh,coeff,tangent)
global OCMATCONT

[failed,infoS,labelS]=OCMATCONT.testadmissibility(tcolmesh,coeff,tangent);

%--< END OF EVALADMISSIBILITYFUNC >--


%---------------------------------
%
% Plot the solution
%
%
function b=PlotContinuationProcess(tcolmesh,coeff,tangent)
global OCMATCONT
b=OCMATCONT.plotcontinuation(tcolmesh,coeff,tangent);

%---------------------------------
%
% Print the solution
%
%
function b=PrintContinuationProcess(tcolmesh,coeff,tangent)
global OCMATCONT
b=OCMATCONT.printcontinuation(tcolmesh,coeff,tangent);

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
% Locate singularity tcolmeshs between x1 and x2
% First locating zeros of testfunctions, then glue them
%
%------------------------------------------------------------

function [tcolmeshs,coeffs,tangents]=LocateSingularity(si)
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
        tcolmeshs=OCMATCONT.testzerosol(idx).tcolmesh;
        coeffs=OCMATCONT.testzerosol(idx).coeff;
        tangents=OCMATCONT.testzerosol(idx).tangent;

    otherwise
        tcolmeshz=zeros(len,OCMATCONT.ndim);
        coeffz=zeros(OCMATCONT.ndim,len);
        tangentz=zeros(OCMATCONT.ndim,len);
        nm=zeros(1,len);

        for ii=1:len
            tcolmeshz(ii,:)=OCMATCONT.testzerosol(idx(ii)).tcolmesh;
            coeffz(:,ii)=OCMATCONT.testzerosol(idx(ii)).coeff;
            tangentz(:,ii)=OCMATCONT.testzerosol(idx(ii)).tangent;
            nm(ii)=norm(tcolmeshz(:,ii));
        end

        if max(nm)-min(nm) < OCMATCONT.OPTIONS.funtolerance
            tcolmeshs=mean(tcolmeshz,1)';
            coeffs=mean(coeffz',1)';
            tangents=mean(tangentz',1)';
        else
            tcolmeshs=[];
            coeffs=[];
            tangents=[];
            return;
        end
end

if lnz==0, return; end

% checking non zeros

DefaultProcessor(tcolmeshs,coeffs,tangents);
tval=EvalTestFunc(nzs,tcolmeshs,coeffs,tangents);
if any(abs(tval(nzs))<=OCMATCONT.OPTIONS.testtolerance)
    tcolmeshs=[];
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

function [tcolmeshs,coeffs,tangents,sing]=xssort(tcolmesh1,tcolmeshs,coeffs,tangents,sing)
% WM: Matlab has a sort function, beter use it...
len=size(tcolmeshs,1);
if len > 1
    tcolmesho=tcolmesh1(ones(1,len),:);
    [dummy,i]=sortrows((tcolmeshs-tcolmesho));
    tcolmeshs=tcolmeshs(i,:);
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

function [failed,f,s]=DefaultProcessor(tcolmesh,coeff,tangent,s)
global OCMATCONT
% WM: this now actually calls the default processor,
% either with or without a singular point structure

if nargin > 3
    [failed,f,s]=OCMATCONT.defaultprocessor(tcolmesh,coeff,tangent,s);
else
    [failed,f]=OCMATCONT.defaultprocessor(tcolmesh,coeff,tangent);
end


%------------------------------------------------
%
% Call method specific mesh adaptation
%
%------------------------------------------------


function [tcolmesh,coeff,tangent,newmesh]=AdaptMesh(tcolmesh,coeff,tangent)
global OCMATCONT
newmesh=false;
[t,y,z,freepar,mdaelpar]=OCMATCONT.drearr(tcolmesh,coeff);
switch OCMATCONT.bvpmethod
    case 'bvp6c'
        res=residual_bvp6c(t,y,freepar,mdaelpar,OCMATCONT.RHS(1:OCMATCONT.HE.numdvariablesmc),OCMATCONT.dae);
        if max(res) < OCMATCONT.OPTIONS.meshadaptreltol;
            return
        end
    case 'bvp5c'
        res=residual_bvp5c(t,y,freepar,mdaelpar,OCMATCONT.RHS(1:OCMATCONT.HE.numdvariablesmc),OCMATCONT.dae);
        if verifyresidual_bvp5c(max(res))
            return
        end
    case 'bvp4c'
        if isempty(OCMATCONT.RHS)
            OCMATCONT.RHS=OCMATCONT.operatoreq(tcolmesh,coeff,[],OCMATCONT.dae,OCMATCONT.bc,OCMATCONT.icfun);
        end
        res=residual_bvp4c(t,y,freepar,mdaelpar,OCMATCONT.RHS(1:OCMATCONT.HE.numdvariablesmc),OCMATCONT.dae);
        if max(res) < OCMATCONT.OPTIONS.meshadaptreltol;
            [tcolmesh,ynew,tangent]=removepoints_bvp4c(t,y,tangent,res);
            coeff=[ynew(:);freepar(:)];
            return
        end
end
% Detect mesh oscillations:  Was there a mesh with
% the same number of ndaes and a similar residual?
% residualReduction = abs(OCMATCONT.meshHistory((OCMATCONT.meshHistory(:,1) == OCMATCONT.HE.TIMEDDATA.nummesh),2) - maxres)/maxres;
oscLikely = false;%any( residualReduction < OCMATCONT.OPTIONS.residualreductionguard);

switch OCMATCONT.bvpmethod
    case 'bvp6c'
        [tcolmesh,ynew,tangent]=meshadaptation_bvp6c(t,y,tangent,res,~oscLikely);
        coeff=[ynew(:);freepar(:)];
    case 'bvp5c'
        [tcolmesh,ynew,tangent]=meshadaptation_bvp5c(t,y,tangent,res,~oscLikely,freepar,mdaelpar,OCMATCONT.dae);
        coeff=[ynew(:)];
    case 'bvp4c'
        [tcolmesh,ynew,tangent]=meshadaptation_bvp4c(t,y,tangent,res,~oscLikely);
        coeff=[ynew(:);freepar(:)];
end
% if ~isempty(OCMATCONT.dataadaptation)
%     OCMATCONT.dataadaptation(tcolmesh,t);
% end
newmesh=true;



function [res,passed]=TestSolutionIntegrity(tcolmesh,coeff)
global OCMATCONT
[t,y,z,freepar,mdaelpar]=OCMATCONT.drearr(tcolmesh,coeff);
switch OCMATCONT.bvpmethod
    case 'bvp6c'
        res=residual_bvp6c(t,y,freepar,mdaelpar,OCMATCONT.RHS(1:OCMATCONT.HE.numdvariablesmc),OCMATCONT.dae);
        passed=max(res) < OCMATCONT.OPTIONS.meshadaptreltol;
    case 'bvp5c'
        res=residual_bvp5c(t,y,freepar,mdaelpar,OCMATCONT.RHS(1:OCMATCONT.HE.numdvariablesmc),OCMATCONT.dae);
        passed=verifyresidual_bvp5c(max(res));
    case 'bvp4c'
        res=residual_bvp4c(t,y,freepar,mdaelpar,OCMATCONT.RHS(1:OCMATCONT.HE.numdvariablesmc),OCMATCONT.dae);
        passed=max(res) < OCMATCONT.OPTIONS.meshadaptreltol;
end

