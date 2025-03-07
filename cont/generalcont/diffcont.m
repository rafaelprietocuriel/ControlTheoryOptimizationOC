function sout=diffcont(varargin)
%
% DIFFCONT main OCMat continuation file for discrete time BVPs
%
% DIFFCONT(PROBLEMTYPE,INITSOL) argument PROBLEMTYPE is a string
% characterizing the type of problem to be solved:
%   'dextremal2fp'   ... saddle-path of an equilibrium, continuing along the
%                       initial point
%   'dextremalp2fp'  ... saddle-path of an equilibrium, continuing along a
%                       parameter value
%   'dextremalt2fp'  ... saddle-path of an equilibrium, continuing the
%                       truncation time
%   'dindifferencesolution'  ... continuation of an indifference threshold
%                       (Skiba point)
%   'limitdextremal' ... continuation of a limitpoint solution
%
% INITSOL is an initial function structure to start the continuation,
% returned by an initialization function, e.g. initocmat_AE_EP,
% initocmat_AE_IS (see OCMat manual)
%
% DIFFCONT(PROBLEMTYPE,INITSOL,INITTANGENT) INITTANGENT is an initial
% tangent, usually this argument is empty.  
%
% DIFFCONT(PROBLEMTYPE,INITSOL,INITTANGENT,OPT) the option structure OPT provides a
% multitude of settings for the continuation process and BVP solver (see
% OCMat manual). 
%
% SOUT=DIFFCONT(...) the structure array consists at least of two elements
% from the intial and last step of the continuation.
%
% The main structure of the diffcont file is taken from the MatCont file 'cont'

global OCMATCONT 
[OCMATCONT.problem_func,sol0,opt]=ParseCommandLine(varargin{:});

delete(findall(0,'Tag','OCMATUserStop')) % delete possible message boxes from a previous continuation

% function handles of the actual continuation type
problemhandles=feval(OCMATCONT.problem_func);
OCMATCONT.operatoreq=problemhandles{1};
OCMATCONT.frechetder=problemhandles{2};
OCMATCONT.findarcposition=problemhandles{6};
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
OCMATCONT.workspaceadapt=problemhandles{19};
OCMATCONT.formatsolution=problemhandles{22};
OCMATCONT.testadmissibility=problemhandles{23};
OCMATCONT.drearr=problemhandles{24};
OCMATCONT.saveintermediate=problemhandles{26};

map=problemhandles{4}{1};
bc=problemhandles{4}{2};
ic=problemhandles{4}{3};
if getocoptions(opt,'SBVPOC','FJacobian')
    mapjac=problemhandles{5}{1};
else
    mapjac=[];
end
if getocoptions(opt,'SBVPOC','BCJacobian')
    bcjac=problemhandles{5}{2};
else
    bcjac=[];
end
if getocoptions(opt,'SBVPOC','ICJacobian')
    icjac=problemhandles{5}{3};
else
    icjac=[];
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

% used options in the continuation file
stepwidth=getocoptions(opt,'OCCONTARG','InitStepWidth'); % initialize step width
stepwidth_max=getocoptions(opt,'OCCONTARG','MaxStepWidth');
stepwidth_min=getocoptions(opt,'OCCONTARG','MinStepWidth');
stepwidth_inc_fac=getocoptions(opt,'OCCONTARG','IncreaseFactor');
stepwidth_dec_fac=getocoptions(opt,'OCCONTARG','DecreaseFactor');
dir_check_step=getocoptions(opt,'OCCONTARG','CheckStep');
dir_check_angle=getocoptions(opt,'OCCONTARG','CheckAngle');

OCMATCONT.OPTIONS.newtonabstol=getocoptions(opt,'NEWTON','AbsTol');
OCMATCONT.OPTIONS.newtonreltol=getocoptions(opt,'NEWTON','RelTol');
OCMATCONT.OPTIONS.maxnewtiter=getocoptions(opt,'NEWTON','MaxNewtonIters');
OCMATCONT.OPTIONS.maxprobes=getocoptions(opt,'NEWTON','MaxProbes');
OCMATCONT.OPTIONS.trm=getocoptions(opt,'NEWTON','TRM');
OCMATCONT.OPTIONS.lambdamin=getocoptions(opt,'NEWTON','LambdaMin');
OCMATCONT.OPTIONS.updatejacfactor=getocoptions(opt,'NEWTON','UpdateJacFactor');
OCMATCONT.OPTIONS.switchtoffnfactor=getocoptions(opt,'NEWTON','SwitchToFFNFactor');
OCMATCONT.OPTIONS.checksingular=logical(strcmp(getocoptions(opt,'NEWTON','CheckSingular'),'on'));
OCMATCONT.OPTIONS.display=getocoptions(opt,'NEWTON','Display');
OCMATCONT.OPTIONS.log=getocoptions(opt,'NEWTON','Log');
OCMATCONT.OPTIONS.admissibletol=getocoptions(opt,'GENERAL','AdmissibleTolerance');
OCMATCONT.newtonsolver=str2func('newtcorr4diffbvp');

OCMATCONT.OPTIONS.xyvectorized=strcmpi(getocoptions(opt,'SBVPOC','Vectorized'),'on');

OCMATCONT.OPTIONS.totalrelativedistance=getocoptions(opt,'OCCONTARG','TotalRelativeDistance');
OCMATCONT.OPTIONS.maxtestiters=getocoptions(opt,'OCCONTARG','MaxTestIters');
OCMATCONT.OPTIONS.testtolerance=getocoptions(opt,'OCCONTARG','TestTolerance');
OCMATCONT.OPTIONS.funtolerance=getocoptions(opt,'OCCONTARG','FunTolerance');
OCMATCONT.OPTIONS.increment=getocoptions(opt,'OCCONTARG','Increment');
OCMATCONT.OPTIONS.vartolerance=getocoptions(opt,'OCCONTARG','VarTolerance');
OCMATCONT.OPTIONS.saveintermediate=strcmp(getocoptions(opt,'OCCONTARG','SaveIntermediate'),'on');
OCMATCONT.OPTIONS.continuationmethod=getocoptions(opt,'OCCONTARG','ContinuationMethod'); % 0 ... Arc-length continuation, 1 Moore-Penrose continuation

% Handle argument functions and additional arguments
[OCMATCONT.map,OCMATCONT.bc,OCMATCONT.ic,OCMATCONT.mapjac,OCMATCONT.bcjac,OCMATCONT.icjac,tmesh0,coeff0,tangent0]= ...
    initialization(sol0,map,bc,ic,mapjac,bcjac,icjac);

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
    OCMATCONT.hitvzero=zeros(2,1);  % tangent where testf is zero
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
    [hitval failed]=EvalHitFunc(1,tmesh0,coeff0,tangent0);
    OCMATCONT.hitvzero(2,:)=hitval;
    if ~isempty(failed)
        error('Evaluation of test functions failed at starting solution.');
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
OCMATCONT.contnum=contnum;

while contnum <= MaxContStep && ~exitflag && ~UserStop.Stop()
    corrections=1;
    meshadaptfailed=0;
    nonadmissible=0;
    while 1 && ~UserStop.Stop()
        % predict next solution
        coeffpre=coeff1 + stepwidth * tangent1;

        [tmesh2,coeff2,tangent2,iter]=OCMATCONT.newtonsolver(tmesh1,coeffpre,tangent1);
        if ~isempty(coeff2) && ((contnum < dir_check_step) || (tangent1'*tangent2 > dir_check_angle))
            nonadmissible=0;
            % test admissibility
            if CheckAdmissibility
                [nonadmissible,infoS,labelS]=EvalAdmissibilityFunc(tmesh2,coeff2,tangent2);
                if nonadmissible
                    % determine new arc structure
                    [out,failed]=NewArcStructure(tmesh2,coeff2,tangent2);
                    ocmatmsg(' \nNon admissible solution detected, reduce stepwidth.\n')
                end
            end
            failed=[];
            if Singularities
                [tfvals,failed]=EvalTestFunc(ActTest,tmesh2,coeff2,tangent2);
                OCMATCONT.testvals(OCMATCONT.atv,:)=tfvals(ActTest);
            end
            if isempty(failed)||~failed
                break
            end
        end
        % reset mesh specific entries of global variable OCMATCONT to previous mesh tmesh1
        if stepwidth > stepwidth_min
            stepwidth=max(stepwidth_min, stepwidth*stepwidth_dec_fac);
            corrections=corrections + 1;
        else      % if not then fail
            if nonadmissible
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
            ocmatmsg('Current step size too small (point %d)\n',contnum);
            exitflag=1;
            if isempty(coeff2)
                tmesh2=tmesh1;
                coeff2=coeff1;
                tangent2=tangent1;
            end
            break
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
    % show information during continuation on command window
    if PrintContStats
        fprintf(1,'\n Continuation step No.: %i\n',contnum);
        fprintf(1,' stepwidth: %g\n',stepwidth);
        fprintf(1,' Newton Iterations: %g\n',iter);
        fprintf(1,' Mesh size: %g\n',length(tmesh2));
        PrintContinuationProcess(tmesh2,coeff2,tangent2);
    end
    if PlotCont
        PlotContinuationProcess(tmesh2,coeff2,tangent2);
    end
    if stepwidth < stepwidth_max && corrections==1 && iter < 4 && ~meshadaptfailed
        stepwidth=min(stepwidth*stepwidth_inc_fac, stepwidth_max);
    end
    if HitTargetValue
        [hitval,failed]=EvalHitFunc(1,tmesh2,coeff2,tangent2);
        if isempty(failed) || ~failed
            OCMATCONT.hitvzero(OCMATCONT.ahv,:)=hitval;
            OCMATCONT.ahv=3-OCMATCONT.ahv;
            testchanges=sign(OCMATCONT.hitvzero(1,:))~=sign(OCMATCONT.hitvzero(2,:));
            if any(testchanges)
                [tmesh2,coeff2,tangent2]=FindSolAtTargetPoint(1,tmesh1,coeff1,tangent1,tmesh2,coeff2,tangent2);
                if ~isempty(coeff2)
                    s.index=contnum;
                    s.label='HTV';
                    s.data.sol=OCMATCONT.formatsolution(tmesh2,coeff2,tangent2);
                    s.msg='This is the solution at the target value.';
                    sout=[sout; s];
                    fprintf('\n Target value hit.\n');
                    fprintf(' label=%s\n Continuation parameter=%g\n', s.label,coeff2(OCMATCONT.HE.contparametercoord));
                    if ExitOnTargetValue
                        exitflag=1;
                    end
                else
                    tmesh2=tmesh1;
                    coeff2=coeff1;
                    tangent2=tangent1;
                    if ExitOnTargetValue
                        exitflag=1;
                    end
                end
            end
        end
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
    if SaveIntermediate
        bvpout(contnum+1).tmesh=tmesh2;
        bvpout(contnum+1).coeff=coeff2;
        bvpout(contnum+1).tangent=tangent2;
        failed=OCMATCONT.saveintermediate(sout,bvpout,contnum+1);
    end
    contnum=contnum+1;
    OCMATCONT.contnum=contnum;
    % shift tmesh1,tangent1  %
    tmesh1=tmesh2;
    coeff1=coeff2;
    tangent1=tangent2;
end

fprintf('\n');
EndTime = clock;

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
global OCMATCONT
coeff=[];
tangent=[];

% no tangent vector given, cycle through base-vectors
ndim=length(coeff0);
ii=ndim;
if ~isempty(tangent0)
    [tmesh,coeff,tangent]=OCMATCONT.newtonsolver(tmesh0,coeff0,tangent0);
end
tangent0=zeros(ndim,1);
while isempty(coeff) && ii>=1
    tangent0(ii)=1;
    try
        [tmesh,coeff,tangent]=OCMATCONT.newtonsolver(tmesh0,coeff0,tangent0);
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


t1(id)=EvalHitFunc(id,tmesh1,coeff1,tangent1);
t2(id)=EvalHitFunc(id,tmesh2,coeff2,tangent2);

ii=1;
tmesh=[];
tangent=[];
coeff=[];
tmax=10*max(abs(t1(id)),abs(t2(id)));
p=1;
try
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
        [tmesh,coeff,tangent]=OCMATCONT.newtonsolver(tmesh3,coeff3,tangent3);

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
catch
    tmesh=[];
    coeff=[];
    tangent=[];
    %ii=OCMATCONT.OPTIONS.maxtestiters;
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
try
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
        [tmesh,coeff,tangent]=OCMATCONT.newtonsolver(tmesh3,coeff3,tangent3);

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
catch
    tmesh=[];
    coeff=[];
    tangent=[];
    ii=OCMATCONT.OPTIONS.maxtestiters;
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
%  Determine new arc structure
%
%------------------------------------------

function [out,failed]=NewArcStructure(tmesh,coeff,tangent)
global OCMATCONT

[out,failed]=feval(OCMATCONT.findarcposition,tmesh,coeff,tangent);

%------------------------------------------
%
%  Evaluate targetvaluefunctions
%
%------------------------------------------

function [out,failed]=EvalHitFunc(id,tmesh,coeff,tangent)
global OCMATCONT
failed=[];
if id==0
    %[out,failed]=feval(OCMATCONT.curve_testf,1:OCMATCONT.nTest, x, tangent);
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
b=OCMATCONT.plotcontinuation(tmesh,coeff,tangent);

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
    sol.solverinfo.tangent=varargin{1};
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

%----------------------------
%
% Initialization of global variables for the boundary value problem
%
%----------------------------

function [mapfinal,bcfinal,icfinal,jacfinal,bcjacfinal,icjacfinal,tmesh,coeff,tangent]=initialization(solinit,map,bc,ic,mapjac,bcjac,icjac)

clear global OCBVP
global OCMATCONT OCBVP

jacfinal=mapjac;
bcjacfinal=bcjac;
icjacfinal=icjac;

% Validate arguments
if ~isfield(solinit,'x')
    msg=sprintf('The field ''x'' not present in SOLINIT.');
    ocmaterror('MATLAB:diffcont:NoXInSolinit','%s',msg);
elseif ~isfield(solinit,'y')
    msg=sprintf('The field ''y'' not present in SOLINIT.');
    ocmaterror('MATLAB:diffcont:NoXInSolinit','%s',msg);
elseif ~isfield(solinit,'parameters')
    msg=sprintf('The field ''parameters'' not present in SOLINIT.');
    ocmaterror('MATLAB:diffcont:NoXInSolinit','%s',msg);
end

if length(solinit.x) < 2
    msg=sprintf('SOLINIT.x must contain at least the two end points.');
    ocmaterror('MATLAB:diffcont:NoXInSolinit','%s',msg);
end

if isempty(solinit.y)
    msg=sprintf('No initial guess provided in SOLINIT.y.');
    ocmaterror('MATLAB:diffcont:NoXInSolinit','%s',msg);
end
if isfield(solinit,'solverinfo') && isfield(solinit.solverinfo,'tangent')
    tangent=solinit.solverinfo.tangent;
else
    tangent=[];
end
if isfield(OCMATCONT,'multipointbvp')
    OCBVP.multipointbvp=OCMATCONT.multipointbvp;
else
    OCBVP.multipointbvp=false;
end
if isfield(OCMATCONT,'sumconstraint')
    OCBVP.sumconstraint=OCMATCONT.sumconstraint;
else
    OCBVP.sumconstraint=0;
end

OCMATCONT.HE.arcindex=arcarg2arcindex(solinit.arcarg);
OCMATCONT.HE.arcarg=solinit.arcarg;

OCMATCONT.HE.numparameter=numel(solinit.parameters);%additional/continuation parameters
OCMATCONT.HE.numparametermc=OCMATCONT.HE.numparameter-OCMATCONT.codimension;%additional parameters minus continuation parameter
OCBVP.npar=OCMATCONT.HE.numparameter; % nmuber of free paramteres
OCBVP.nparmc=OCMATCONT.HE.numparameter-OCMATCONT.codimension; % number of free paramteres excluded the number of continuation parameters
nummap=size(solinit.y,1);
tmesh=solinit.x;
coeff=solinit.y(:);
coeff=[coeff;solinit.parameters(:)];
OCMATCONT.HE.numdvariables=numel(coeff);
OCMATCONT.HE.numdvariablesmc=OCMATCONT.HE.numdvariables-OCMATCONT.codimension;

OCBVP.nummap=nummap;
OCBVP.N=numel(tmesh);
OCBVP.Np1=OCBVP.N+1;
OCBVP.neqn=OCBVP.nummap;
OCBVP.explicitparameterdependence=true;
OCBVP.nN=OCBVP.nummap*OCBVP.N;
OCBVP.cols=1:OCBVP.nummap;             % in the global Jacobian

OCMATCONT.HE.DDATA.meshvalcoord=reshape(1:OCBVP.N*OCBVP.neqn,OCBVP.neqn,OCBVP.N);
OCMATCONT.HE.parametercoord=numel(OCMATCONT.HE.DDATA.meshvalcoord)+(1:OCMATCONT.HE.numparameter).';
OCMATCONT.HE.contparametercoord=numel(OCMATCONT.HE.DDATA.meshvalcoord)+OCMATCONT.HE.numparameter-1+(1:OCMATCONT.codimension);

OCBVP.numarc=OCMATCONT.numarc;

OCBVP.multiarccalc=false;

if ~OCBVP.multipointbvp
    OCBVP.nBCs=nummap+OCMATCONT.HE.numparameter-OCMATCONT.codimension;
else
    OCBVP.nBCs=nummap*OCMATCONT.multipointorder+OCMATCONT.HE.numparameter-OCMATCONT.codimension-OCMATCONT.numsumconstraint;
    OCBVP.nICs=OCMATCONT.numsumconstraint;
end
OCBVP.rows=OCBVP.nBCs+1:OCBVP.nBCs+OCBVP.nummap;   % define the action area
OCBVP.Lidx=solinit.arcposition(1,:);
OCBVP.Ridx=solinit.arcposition(2,:);
OCBVP.Nint=OCBVP.Ridx-OCBVP.Lidx;
OCBVP.arcposition=solinit.arcposition;
if isfield(OCMATCONT,'multipointbvp') && OCMATCONT.multipointbvp
    OCBVP.pathcoord=solinit.pathposition;
    OCBVP.multiarccalc=true;
    OCBVP.numpath=length(OCBVP.pathcoord(1,:));
else
    OCBVP.pathcoord=[OCBVP.Lidx(1);OCBVP.Ridx(end)];
    OCBVP.multiarccalc=false;
    OCBVP.numpath=1;
end
% Function handles
mapfinal=@(x,y,region,p,par)map(x,y,region,p,par);
if isa(mapjac,'function_handle')
    jacfinal=@(x,y,region,p,par)mapjac(x,y,region,p,par);
end
bcfinal=@(ya,yb,p,par) bc(ya,yb,p,par);
if isa(bcjac,'function_handle')
    bcjacfinal=@(ya,yb,p,par) bcjac(ya,yb,p,par);
end
if isfield(OCMATCONT,'sumconstraint')
    icfinal=@(x,y,p,par) ic(x,y,p,par);
    if isa(icjac,'function_handle')
        icjacfinal=@(x,y,p,par) icjac(x,y,p,par);
    end
else
    icfinal=[];
    icjacfinal=[];
end
% options for numerical differentiation
Joptions=[];
dPoptions=[];
dBCoptions=[];
dICoptions=[];
dICdPoptions=[];
if 1%isempty(jacfinal)
    Joptions.diffvar=2;  % dF(x,y)/dy
    if OCMATCONT.OPTIONS.xyvectorized
        Joptions.vectvars=[1,2];
    else
        Joptions.vectvars=[];
    end
    dPoptions.diffvar=6;  % dF(x,y,region,p)/dp
    dPoptions.vectvars=[]; % no vectorization for parameter Jacobian
end
if isempty(bcjacfinal)
    dBCoptions.vectvars=[];
end
if 1%isempty(icjacfinal)
    dICoptions.diffvar=2;
    if OCMATCONT.OPTIONS.xyvectorized
        %dICoptions.vectvars=[1,2];
        dICoptions.vectvars=[];
    else
        dICoptions.vectvars=[];
    end
    dICdPoptions.diffvar=5;  % dF(x,y,region,p)/dp
    dICdPoptions.vectvars=[]; % no vectorization for parameter Jacobian
end
OCBVP.Joptions=Joptions;
OCBVP.dPoptions=dPoptions;
OCBVP.dBCoptions=dBCoptions;
OCBVP.dICoptions=dICoptions;
OCBVP.dICdPoptions=dICdPoptions;
OCBVP.averageJac=isempty(jacfinal);

warnoffId={'MATLAB:singularMatrix','MATLAB:nearlySingularMatrix'};
for ii=1:length(warnoffId)
    warnstat(ii)=warning('query',warnoffId{ii});
    warnoff(ii)=warnstat(ii);
    warnoff(ii).state='off';
end
OCBVP.warnstat=warnstat;
OCBVP.warnoff=warnoff;
OCBVP.warnoffId=warnoffId;
