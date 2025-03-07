function [sout bvpout]=bvpconts(varargin)

global OCBVP OCMATCONT

[OCMATCONT.problem_func,solinit,opt]=ParseCommandLine(varargin{:});

delete(findall(0,'Tag','OCMATUserStop')) % delete possible message boxes from a previous continuation

%
OCMATCONT.bvpmethod=getocoptions(opt,'GENERAL','BVPMethod');
bvpsolver=str2func(OCMATCONT.bvpmethod);

OCMATCONT.OPTIONS.totalrelativedistance=getocoptions(opt,'OCCONTARG','TotalRelativeDistance');
OCMATCONT.OPTIONS.maxtestiters=getocoptions(opt,'OCCONTARG','MaxTestIters');
OCMATCONT.OPTIONS.testtolerance=getocoptions(opt,'OCCONTARG','TestTolerance');
OCMATCONT.OPTIONS.funtolerance=getocoptions(opt,'OCCONTARG','FunTolerance');
OCMATCONT.OPTIONS.increment=getocoptions(opt,'OCCONTARG','Increment');
OCMATCONT.OPTIONS.vartolerance=getocoptions(opt,'OCCONTARG','VarTolerance');

SaveIntermediate=strcmpi(getocoptions(opt,'OCCONTARG','SaveIntermediate'),'on');
MaxContStep=getocoptions(opt,'OCCONTARG','MaxContinuationSteps');
PrintContStats=strcmpi(getocoptions(opt,'OCCONTARG','PrintContStats'),'on');
PlotCont=strcmpi(getocoptions(opt,'OCCONTARG','PlotCont'),'on');
HitTargetValue=getocoptions(opt,'OCCONTARG','HitTargetValue');
ExitOnTargetValue=strcmpi(getocoptions(opt,'OCCONTARG','ExitOnTargetValue'),'on');

stepwidth=getocoptions(opt,'OCCONTARG','InitStepWidth'); % initialize step width
stepwidth_max=getocoptions(opt,'OCCONTARG','MaxStepWidth');
stepwidth_min=getocoptions(opt,'OCCONTARG','MinStepWidth');
stepwidth_inc_fac=getocoptions(opt,'OCCONTARG','IncreaseFactor');
stepwidth_dec_fac=getocoptions(opt,'OCCONTARG','DecreaseFactor');

problemhandles=feval(OCMATCONT.problem_func);
OCMATCONT.targetvaluefunc=problemhandles{9};
OCMATCONT.plotcontinuation=problemhandles{21};
OCMATCONT.formatsolution=problemhandles{22};
OCMATCONT.saveintermediate=problemhandles{26};
OCMATCONT.printcontinuation=problemhandles{30};
if getocoptions(opt,'SBVPOC','FJacobian')
    odejac=problemhandles{5}{1};
else
    odejac=[];
end
if getocoptions(opt,'SBVPOC','BCJacobian')
    bcjac=problemhandles{5}{2};
else
    bcjac=[];
end

[solinit,freepar,odefinal,bcfinal,jacfinal,bcjacfinal]=adapt2solver(OCMATCONT.bvpmethod,solinit,problemhandles{4}{1},problemhandles{4}{2},odejac,bcjac,size(solinit.y,1),length(solinit.parameters));


optbvp=opt.BVP;
optbvp=bvpset(optbvp,'RelTol',getocoptions(opt,'SBVPOC','MeshAdaptRelTol'));
optbvp=bvpset(optbvp,'AbsTol',getocoptions(opt,'SBVPOC','MeshAdaptAbsTol'));
optbvp=bvpset(optbvp,'Vectorized',getocoptions(opt,'SBVPOC','Vectorized'));
optbvp=bvpset(optbvp,'FJacobian',jacfinal);
optbvp=bvpset(optbvp,'BCJacobian',bcjacfinal);


if HitTargetValue
    OCMATCONT.hitvzero=zeros(2,1);  % tangent where testf is zero
    OCMATCONT.ahv=1;
end

if PlotCont
    clf
end
StartTime = clock;

contnum=0;

UserStop=stoploop('Stop continuation.');
sout=[];
bvpout=[];
solold=[];

while contnum<=MaxContStep && ~UserStop.Stop()
    contnum=contnum+1;
    ti1=tic;
    solinit=transform2matlabform(OCMATCONT.bvpmethod,solinit);
    try
        [sol solinfo]=bvpsolver(odefinal,bcfinal,solinit,optbvp);
    catch ME
        fprintf(1,'%s\n',ME.identifier)
        sol.err=1;
    end
    sol=transform2ocmatform(OCMATCONT.bvpmethod,sol);
    ti1=toc(ti1);
    if contnum==1
        if sol.err
            ocmatmsg('No convergence at starting function sol.\n')
            ocmatmsg('elapsed time=%.1f secs\n', etime(clock, StartTime));
            ocmatmsg('0 npoints\n');
            return;
        end
        s.index=1;
        s.label='00';
        s.data.sol=OCMATCONT.formatsolution(sol.x,sol.y(:),OCMATCONT.InitialSecant);
        s.msg='This is the first solution of the BVP continuation';

        fprintf('first solution found\n');
        sout=s;
    end
    alpha=1;
    if sol.err
        contnum=contnum-1;
        stepwidth=stepwidth*stepwidth_dec_fac;
        sol=solold;
        if stepwidth<stepwidth_min
            break
        end
        if ~isempty(solold2)
            solinit=sol;
            OCMATCONT.InitialSecant=sol.y(OCMATCONT.LastIntialDistributionIndex)-solold2.y(OCMATCONT.LastIntialDistributionIndex);
            OCMATCONT.InitialSecantO=OCMATCONT.InitialSecant;
            absv=abs(OCMATCONT.InitialSecant)./max(1e-5,abs(sol.y(OCMATCONT.LastIntialDistributionIndex)));
            endval=OCMATCONT.InitialSecant(end);
            OCMATCONT.InitialSecant(max(absv)/2>absv)=0;
            OCMATCONT.InitialSecant(end)=endval;
            numpar=length(OCMATCONT.HE.parametercoord);
            OCMATCONT.InitialSecant(end-numpar+1:end-1)=0;

            OCMATCONT.InitialSecant=stepwidth*OCMATCONT.InitialSecant/norm(OCMATCONT.InitialSecant);
            OCMATCONT.LastIntialDistribution=sol.y(OCMATCONT.LastIntialDistributionIndex);
        else
            solinit=sol;
            OCMATCONT.LastIntialDistribution=solinit.y(OCMATCONT.LastIntialDistributionIndex);
            OCMATCONT.InitialSecant=[zeros(OCBVP.numode+OCBVP.numparameter-1,1);stepwidth];
        end
    else
        if solinfo.itnl<4
            oldstepwidth=stepwidth;
            stepwidth=stepwidth*stepwidth_inc_fac;
            if stepwidth>stepwidth_max
                stepwidth=stepwidth_max;
            end
            %alpha=alpha*stepwidth/oldstepwidth;
        end
        if ~isempty(solold)
            if length(solold.x)==length(sol.x) && all(solold.x-sol.x)==0
                solinit.x=sol.x;
                solinit.y=sol.y+alpha*(sol.y-solold.y);
            else
                solold=evalatsamemesh(OCMATCONT.bvpmethod,solold,sol);
                solinit.x=sol.x;
                solinit.y=sol.y+alpha*(sol.y-solold.y);
            end
            OCMATCONT.InitialSecant=sol.y(OCMATCONT.LastIntialDistributionIndex)-solold.y(OCMATCONT.LastIntialDistributionIndex);
            OCMATCONT.InitialSecantO=OCMATCONT.InitialSecant;
                       absv=abs(OCMATCONT.InitialSecant)./max(1e-5,abs(sol.y(OCMATCONT.LastIntialDistributionIndex)));
                        endval=OCMATCONT.InitialSecant(end);
                        OCMATCONT.InitialSecant(max(absv)/2>absv)=0;
                        OCMATCONT.InitialSecant(end)=endval;
            numpar=length(OCMATCONT.HE.parametercoord);
            OCMATCONT.InitialSecant(end-numpar+1:end-1)=0;
            OCMATCONT.InitialSecant=stepwidth*OCMATCONT.InitialSecant/norm(OCMATCONT.InitialSecant);
            OCMATCONT.LastIntialDistribution=sol.y(OCMATCONT.LastIntialDistributionIndex);
        else
            solinit=sol;
            OCMATCONT.LastIntialDistribution=solinit.y(OCMATCONT.LastIntialDistributionIndex);
            OCMATCONT.InitialSecant=[zeros(OCBVP.numode+OCBVP.numparameter-1,1);stepwidth];
        end

        if PrintContStats
            fprintf(1,'\n Continuation step No.: %i\n',contnum);
            fprintf(1,' stepwidth: %g\n',stepwidth);
            if length(solinfo.itlin)==2
                fprintf(1,' Calculation Time: %4.3f\n Newton Iterations: linear: [%d, %d],\tnonlinear: %d\n',ti1,solinfo.itlin(1),solinfo.itlin(2),solinfo.itnl)
            elseif length(solinfo.itlin)==2
                fprintf(1,' Calculation Time: %4.3f\n Newton Iterations: linear: %d,\tnonlinear: %d\n',ti1,solinfo.itlin(1),solinfo.itnl)
            else
                fprintf(1,' Calculation Time: %4.3f\n Newton Iterations: linear: %d,\tnonlinear: %d\n',ti1,0,solinfo.itnl)
            end
            fprintf(1,' Mesh size: %g\n',length(sol.x));
            PrintContinuationProcess(sol.x,sol.y(:),[]);
        end
        if PlotCont
            OCMATCONT.plotcontinuation(sol.x,sol.y(:));
        end
        if HitTargetValue
            [hitval,failed]=EvalHitFunc(1,sol.x,sol.y(:),[]);
            if isempty(failed)|| ~failed
                OCMATCONT.hitvzero(OCMATCONT.ahv,:)=hitval;
                OCMATCONT.ahv=3-OCMATCONT.ahv;
                testchanges=sign(OCMATCONT.hitvzero(1,:))~=sign(OCMATCONT.hitvzero(2,:));
                if any(testchanges) && contnum>1
                    [solh solinfoh]=FindSolAtTargetPoint(1,solold,sol);
                    if ~isempty(solh)
                        contnum=contnum+1;
                        s.index=contnum;
                        s.label='HTV';
                        s.data.sol=OCMATCONT.formatsolution(solh.x,solh.y,[]);
                        s.msg='This is the solution at the target value.';
                        sout=[sout; s];
                        fprintf('\n Target value hit.\n');
                        fprintf(' label=%s\n Continuation parameter=%g\n', s.label,s.data.sol.continuationparameter);
                        if SaveIntermediate
                            bvpout(contnum).tmesh=solh.x;
                            bvpout(contnum).coeff=solh.y(:);
                            bvpout(contnum).secant=OCMATCONT.InitialSecant;
                            bvpout(contnum).bvpinfo=solinfoh;
                            failed=OCMATCONT.saveintermediate(sout,bvpout,1);
                        end
                        if ExitOnTargetValue
                            break
                        end
                    else
                        if ExitOnTargetValue
                            break
                        end
                    end
                end
            end
        end

        solold2=solold;
        solold=sol;
        if SaveIntermediate
            bvpout(contnum).tmesh=sol.x;
            bvpout(contnum).coeff=sol.y(:);
            bvpout(contnum).secant=OCMATCONT.InitialSecant;
            bvpout(contnum).bvpinfo=solinfo;
            failed=OCMATCONT.saveintermediate(sout,bvpout,1);
        end
    end
end
fprintf('\n');
EndTime = clock;

s.index=contnum;
s.label='99';
s.data.sol=OCMATCONT.formatsolution(sol.x,sol.y(:),OCMATCONT.InitialSecant);
s.msg='This is the last solution of the BVP continuation';
sout=[sout; s];
if SaveIntermediate
    bvpout(contnum+1).tmesh=sol.x;
    bvpout(contnum+1).coeff=sol.y;
    bvpout(contnum+1).tangent=OCMATCONT.InitialSecant;
    bvpout(contnum+1).bvpinfo=solinfo;
    failed=OCMATCONT.saveintermediate(sout,bvpout,contnum+1);
end

UserStop.Clear();
clear UserStop
fprintf('elapsed time  = %.1f secs\n', etime(EndTime, StartTime));

%-------------------------------------
%
% Find the solution at the target point
%
%-------------------------------------
    function [sol solinfo]=FindSolAtTargetPoint(id,sol1,sol2)
        OCMATCONT.LastIntialDistribution0=OCMATCONT.LastIntialDistribution;
        OCMATCONT.InitialSecant0=OCMATCONT.InitialSecant;
        t1(id)=EvalHitFunc(id,sol1.x,sol1.y(:),[]);
        t2(id)=EvalHitFunc(id,sol2.x,sol2.y(:),[]);

        ii=1;
        tmax=10*max(abs(t1(id)),abs(t2(id)));
        try
            while ii<=OCMATCONT.OPTIONS.maxtestiters
                OCMATCONT.LastIntialDistribution=sol1.y(OCMATCONT.LastIntialDistributionIndex,1);
                OCMATCONT.InitialSecant=[zeros(OCBVP.numode,1);t1];
                sol=transform2matlabform(OCMATCONT.bvpmethod,sol1);
                try
                    [sol solinfo]=bvpsolver(odefinal,bcfinal,sol,optbvp);
                catch ME
                    fprintf(1,'%s\n',ME.identifier)
                    sol.err=1;
                end
                sol=transform2ocmatform(OCMATCONT.bvpmethod,sol);
                if sol.err
                    OCMATCONT.LastIntialDistribution=sol2.y(OCMATCONT.LastIntialDistributionIndex,1);
                    OCMATCONT.InitialSecant=[zeros(OCBVP.numode,1);t2];
                    sol=transform2matlabform(OCMATCONT.bvpmethod,sol2);
                    try
                        [sol solinfo]=bvpsolver(odefinal,bcfinal,sol,optbvp);
                    catch ME
                        fprintf(1,'%s\n',ME.identifier)
                        sol.err=1;
                    end
                    sol=transform2ocmatform(OCMATCONT.bvpmethod,sol);
                end

                if sol.err
                    return
                end
                tval=EvalHitFunc(id,sol.x,sol.y(:),[]);
                if abs(tval(id))>tmax
                    sol=[];
                    OCMATCONT.LastIntialDistribution=OCMATCONT.LastIntialDistribution0;
                    OCMATCONT.InitialSecant=OCMATCONT.InitialSecant0;
                    break;
                end
                if abs(tval(id))<=OCMATCONT.OPTIONS.testtolerance
                    OCMATCONT.LastIntialDistribution=OCMATCONT.LastIntialDistribution0;
                    OCMATCONT.InitialSecant=OCMATCONT.InitialSecant0;
                    break;
                end
                ii=ii+1;
                sol=[];
            end
        catch
            OCMATCONT.LastIntialDistribution=OCMATCONT.LastIntialDistribution0;
            OCMATCONT.InitialSecant=OCMATCONT.InitialSecant0;
            sol=[];
        end
    end
end

function [solinit,freepar,odefinal,bcfinal,jacfinal,bcjacfinal]=adapt2solver(bvpmethod,solinit,odefun,bcfun,jac,bcjac,neqn,nparam)

global OCBVP OCMATCONT


% make functions feval-free
ode = fcnchk(odefun);
bc  = fcnchk(bcfun);
if ischar(jac)
    jac = str2func(jac);
end
if ischar(bcjac)
    bcjac = str2func(bcjac);
end

jacfinal = jac;
bcjacfinal = bcjac;

switch bvpmethod
    case {'mtom0','tom'}
        numparam=length(solinit.parameters);
        numode=size(solinit.y,1);
        OCBVP.numparameter=numparam;
        OCBVP.numode=numode;
        N=length(solinit.x);

        OCMATCONT.HE.DDATA.meshvalcoord=reshape(1:N*(numode+numparam),numode+numparam,N);
        OCMATCONT.HE.DDATA.meshvalcoord(end-numparam+1:end,:)=[];
        OCMATCONT.HE.parametercoord=numode+(1:numparam).';
        OCMATCONT.HE.contparametercoord=OCMATCONT.HE.parametercoord(end);%-1+(1:OCMATCONT.codimension);
        OCMATCONT.LastIntialDistributionIndex=[OCMATCONT.HE.DDATA.meshvalcoord(:,1);OCMATCONT.HE.parametercoord(:)];
        freepar=solinit.parameters;
        OCMATCONT.HE.arcarg=solinit.arcarg;
        OCMATCONT.HE.numarc=length(OCMATCONT.HE.arcarg);
        ff=find(diff(solinit.x)==0);
        OCMATCONT.HE.TIMEDDATA.leftarcindex=[1 ff];
        OCMATCONT.HE.TIMEDDATA.rightarcindex=[ff+1 length(solinit.x)];
            
        solinit.y=[solinit.y;solinit.parameters(ones(1,size(solinit.y,2)),:).'];
        % incorporate unknown parameters (for TOM)
        if ~isfield(solinit.solverinfo,'tangent') || isempty(solinit.solverinfo.tangent)
        OCMATCONT.InitialSecant=[zeros(OCBVP.numode,1);zeros(numparam-1,1);1];
        else
        OCMATCONT.InitialSecant=solinit.solverinfo.tangent(:);
        end
        OCMATCONT.LastIntialDistribution=solinit.y(:,1);
        
        odefinal = @odeParameters;
        if ~isempty(jac)
            jacfinal = @jacParameters;
        end
        bcfinal  = @bcParameters;
        if ~isempty(bcjac)
            bcjacfinal = @bcjacParameters;
        end
    case {'bvp4c','bvp6c'}
        numparam=length(solinit.parameters);
        numode=size(solinit.y,1);
        OCBVP.numparameter=numparam;
        OCBVP.numode=numode;
        N=length(solinit.x);

        OCMATCONT.HE.DDATA.meshvalcoord=reshape(1:N*(numode),numode,N);
        OCMATCONT.HE.parametercoord=numel(OCMATCONT.HE.DDATA.meshvalcoord)+(1:numparam).';
        OCMATCONT.HE.contparametercoord=OCMATCONT.HE.parametercoord(end);%-1+(1:OCMATCONT.codimension);
        OCMATCONT.LastIntialDistributionIndex=[OCMATCONT.HE.DDATA.meshvalcoord(:,1);OCMATCONT.HE.parametercoord(:)];

        freepar=solinit.parameters;
        OCMATCONT.HE.arcarg=solinit.arcarg;
        OCMATCONT.HE.numarc=length(OCMATCONT.HE.arcarg);
        ff=find(diff(solinit.x)==0);
        OCMATCONT.HE.TIMEDDATA.leftarcindex=[1 ff];
        OCMATCONT.HE.TIMEDDATA.rightarcindex=[ff+1 length(solinit.x)];
        
        if ~isfield(solinit.solverinfo,'tangent') || isempty(solinit.solverinfo.tangent)
            OCMATCONT.InitialSecant=[zeros(OCBVP.numode,1);zeros(numparam-1,1);1];
        else
            OCMATCONT.InitialSecant=solinit.solverinfo.tangent(:);
        end
        OCMATCONT.LastIntialDistribution=[solinit.y(:,1);solinit.parameters(:)];
        solinit.y=[solinit.y(:);solinit.parameters(:)];
        odefinal=ode;
        bcfinal=bc;
    otherwise
end

% ---------------------------------------------------------
% Nested functions
% ---------------------------------------------------------

    function f = odeParameters(x,y)
        p = y(neqn+1:neqn+nparam,1);  % extract p from y(:,1)
        % add trivial equations for unknown parameters
        f = [ ode(x,y(1:neqn,:),p);
            zeros(nparam,numel(x))];
    end  % odeParameters

% ---------------------------------------------------------

    function res = bcParameters(ya,yb)
        p = ya(neqn+1:neqn+nparam);   % extract p from ya
        res = bc(ya(1:neqn),yb(1:neqn),p);
    end  % bcParameters

% ---------------------------------------------------------

    function J = jacParameters(x,y)
        p = y(neqn+1:neqn+nparam,1);   % extract p from y(:,1)
        [dfdy,dfdp] = jac(x,y(1:neqn,:),p);
        % add trivial equations for unknown parameters
        J = [dfdy, dfdp;
            zeros(nparam,neqn+nparam)];
    end  % jacParameters
end

%---------------------------------
%
% Print the solution
%
%
function b=PrintContinuationProcess(tmesh,coeff,tangent)
global OCMATCONT
b=OCMATCONT.printcontinuation(tmesh,coeff,tangent);
end


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
end

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
end

function sol=transform2ocmatform(bvpmethod,sol)

switch bvpmethod
    case {'mtom0','tom'}
        sol.y=sol.y(:);
        sol.x=sol.x/sol.x(end);
    case {'bvp4c','bvp6c'}
       sol.y=[sol.y(:);sol.parameters(:)];
end
end

function sol=transform2matlabform(bvpmethod,sol)
global OCMATCONT OCBVP

switch bvpmethod
    case {'mtom0','tom'}
        sol.y=reshape(sol.y,OCBVP.numode+OCBVP.numparameter,[]);
    case {'bvp4c','bvp6c'}
        sol.parameters=sol.y(OCMATCONT.HE.parametercoord);
        sol.y=sol.y(OCMATCONT.HE.DDATA.meshvalcoord);
end
end

function solold=evalatsamemesh(bvpmethod,solold,sol)
global OCMATCONT OCBVP
solold=transform2matlabform(bvpmethod,solold);
OCBVP.N=length(sol.x);
OCBVP.nN=OCBVP.numode*OCBVP.N;
switch bvpmethod
    case {'mtom0','tom'}
        solold.y=interp1(solold.x,solold.y.',sol.x).';
        solold.y=solold.y(:);
        OCMATCONT.HE.DDATA.meshvalcoord=reshape(1:OCBVP.N*(OCBVP.numode+OCBVP.numparameter),OCBVP.numode+OCBVP.numparameter,OCBVP.N);
        OCMATCONT.HE.DDATA.meshvalcoord(end,:)=[];

    case {'bvp4c','bvp6c'}
        [solold.y solold.yp]=deval6c(solold,sol.x);
        solold.x=sol.x;
        solold=transform2ocmatform(bvpmethod,solold);
        OCMATCONT.HE.DDATA.meshvalcoord=reshape(1:OCBVP.N*(OCBVP.numode),OCBVP.numode,OCBVP.N);
        OCMATCONT.HE.parametercoord=numel(OCMATCONT.HE.DDATA.meshvalcoord)+(1:OCBVP.numparameter).';
        OCMATCONT.HE.contparametercoord=OCMATCONT.HE.parametercoord(end);%-1+(1:OCMATCONT.codimension);
        OCMATCONT.LastIntialDistributionIndex=[OCMATCONT.HE.DDATA.meshvalcoord(:,1);OCMATCONT.HE.parametercoord(:)];
        ff=find(diff(sol.x)==0);
        OCMATCONT.HE.TIMEDDATA.leftarcindex=[1 ff];
        OCMATCONT.HE.TIMEDDATA.rightarcindex=[ff+1 OCBVP.N];
end

end