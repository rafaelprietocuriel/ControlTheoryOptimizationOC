function out = finiteimpulse
%
% Equilibrium curve definition file for a problem in odefile
%

out{1}  = @curve_funcnew;
out{2}  = @defaultprocessor;
out{3}  = @options;
out{4}  = @jacobian;
out{5}  = @hessians;
out{6}  = @testf;
out{7}  = @userf;
out{8}  = @process;
out{9}  = @singmat;
out{10} = @locate;
out{11} = @init;
out{12} = @done;
out{13} = [];%@adapt;
out{21}=@saveintermediate;

%-------------------------------------------------------
function res = curve_funcnew(arg)
global ics
arcarg=ics.arcarg;
arcnum=ics.arcnum;
jumparg=ics.jumparg;

%%%%%
[x,p,T,initialstates] = rearrnew(arg);
dynVar=x(ics.genDepVarCoordinates);
timepoints=[ics.InitialTime x(ics.JumpTimeCoordinates).' T];

res=[];
for ii=1:arcnum+1
    dynVarLR=dynVar(:,(2*ii-1):2*ii); % variable representing the left and right limits
    if ii==1
        % initial impulse or continuity
        res=ics.jump(timepoints(ii),dynVarLR,p,jumparg(ii));
        % initial condition
        res=[res; ...
            ics.initialcond(timepoints(ii),dynVarLR,p,arcarg(1),initialstates)];
    elseif ii==arcnum+1
        % end impulse or continuity
        res=[res; ...
            ics.jump(timepoints(ii),dynVarLR,p,jumparg(ii))];
        % transversality condition
        res=[res;ics.transversality(timepoints(ii),dynVarLR,p,arcarg(end))];
    end
    if arcnum>1 && ii>=2 && ii<=arcnum
        % interior jump or continuity
        res=[res; ...
            ics.jump(timepoints(ii),dynVarLR,p,jumparg(ii))];
        if jumparg(ii)
            % optimality of interior impulse
            res=[res; ...
                ics.interiorjumpcond(timepoints(ii),dynVarLR,p,jumparg(ii),arcarg(ii-1))];
        end
    end
    if ii<=arcnum
        %Gamma_i
        % explicit path between impulses
        dynVarI=dynVar(:,2*ii:(2*ii+1)); % variable representing the points on the interval between two impulses
        res=[res;ics.continousdyn(timepoints(ii:ii+1),dynVarI,p,arcarg(ii))];
    end
end
%---------------------------------------------------------------
function jac = jacobian(varargin)
global eds
[x,p] = rearr(varargin{1});
p = num2cell(p);
%varargout{1}=DS_jac('DS_jac',x,p,T);
%jac = [cjac(eds.func,eds.Jacobian,x,p,eds.ActiveParams) cjacp(eds.func,eds.JacobianP,x,p,eds.ActiveParams)];
%jac = [cjac(@thequations,eds.Jacobian,x,p,eds.ActiveParams) cjacp(@thequations,eds.JacobianP,x,p,eds.ActiveParams)];
jac = [cjac(@curve_funcnew,eds.Jacobian,x) cjacp(@curve_funcnew,eds.JacobianP,x)];
jac=sparse(jac);
%---------------------------------------------------------------
function hess = hessians(varargin)
global cds eds
xo = varargin{1}; [x,p] =  rearr(xo);p=num2cell(p);
hh = chess(eds.func,eds.Jacobian,eds.Hessians,x,p,eds.ActiveParams);
hp = chessp(eds.func,eds.Jacobian,eds.HessiansP,x,p,eds.ActiveParams);
x1 = xo; x1(cds.ndim) = x1(cds.ndim) - option.Increment;
x2 = xo; x2(cds.ndim) = x2(cds.ndim) + option.Increment;
hpp = (cjac(cds.curve_func,cds.curve_jacobian,x2,[]) - cjac(cds.curve_func,cds.curve_jacobian,x1,[])) / (2*option.Increment);
for i = 1:cds.ndim-1
    hess(:,:,i) = [ hh(:,:,i) hpp(:,i)];
end
hess(:,:,cds.ndim) = [ hp(:,:) hpp(:,cds.ndim)];
%---------------------------------------------------------------
function varargout = defaultprocessor(varargin)
global eds cds ics
arcarg=ics.arcarg;
jumparg=ics.jumparg;

%%%%%
[x,par,T] = rearrnew(varargin{1});
if nargin > 2
    s = varargin{3};
    s.data.v = eds.v;
    varargout{3} = s;
end
% n=10;
% depvar=x(ics.genDepVarCoordinates);
% timepoints=[ics.InitialTime x(ics.JumpTimeCoordinates).' T];
% out=[];
% for ii=1:ics.arcnum
%     t=linspace(timepoints(ii),timepoints(ii+1),n);
%     tmp=ics.impulseconstraint(t,depvar(:,2*ii),par,arcarg(ii));
%     if ii<=ics.arcnum
%         tmp([1 end])=[];
%     end
% %     if ii==ics.arcnum
% %         fprintf('time : %f\tjumpcd : %f\n',timepoints(ii+1),tmp(end))
% %     end
%     out=[out tmp];
% end
% check integrity
varargout{2}=NaN;%out(:);
% all done succesfully
varargout{1} = 0;%any(abs(diff(sign(out)))>0);
%-------------------------------------------------------------
function option = options
global eds cds
option = contset;
if cds.ndim<3
    option=contset(option,'IgnoreSingularity',2);
end
% Check for symbolic derivatives in odefile

symord = 0;
if ~isempty(eds.Jacobian), symord = 1; end
if ~isempty(eds.Hessians), symord = 2; end
if ~isempty(eds.Der3), symord = 3; end
if ~isempty(eds.Der4), symord = 4; end
if ~isempty(eds.Der5), symord = 5; end

option = contset(option, 'SymDerivative', symord);
option = contset(option, 'Workspace', 0);
%option = contset(option, 'Singularities', 0);
option = contset(option, 'Adapt', 0);

UserInfo(1).name='reachtimehorizon';UserInfo(1).state=1;UserInfo(1).label='IOC_TH'; % becomes zero if the specified time horizon (GlobalImpulseVariable.TimeHorizon) is reached during continuation
UserInfo(2).name='jumpingtimesvstimehorizon';UserInfo(2).state=1;UserInfo(2).label='IOC_JLH'; % becomes zero if a jumping time exceeds actual time horizon
UserInfo(3).name='negativetime';UserInfo(3).state=1;UserInfo(3).label='IOC_NegT'; % the jumping times and timehorizon are monitored not to become negative
UserInfo(4).name='jumpfuncAt_T';UserInfo(4).state=1;UserInfo(4).label='IOC_JCT'; %
UserInfo(5).name='reachparameter';UserInfo(5).state=1;UserInfo(5).label='IOC_RP';
UserInfo(6).name='jumpfuncAt_t0';UserInfo(6).state=1;UserInfo(6).label='IOC_JC0';
%UserInfo(7).name='interiorjumpfunc';UserInfo(7).state=1;UserInfo(7).label='IOC_IJC'; % becomes zero if interior jump condition is satisfied

option = contset(option, 'UserfunctionsInfo',UserInfo);
option = contset(option, 'Userfunctions',1);
%option = contset(option, 'Userfunctions',0);
%option = contset(option,'UserfunctionsStop',[1 1 0 0 0 0 0]);
option = contset(option,'UserfunctionsStop',[]);

symordp = 0;
if ~isempty(eds.JacobianP), symordp = 1; end
if ~isempty(eds.HessiansP),  symordp = 2; end
option = contset(option, 'SymDerivativeP', symordp);

%-----------------------------------------------------------------
function [out, failed] = userf(userinf, id, x, v)
global eds
dim =size(id,2);
failed = [];
out(dim) = 0;
for i=1:dim
    lastwarn('');
    [genDynVar,pararg,T,initialstates] = rearrnew(x);
    pararg = num2cell([initialstates pararg T]);
    if (userinf(i).state==1)
        out(i)=feval(eds.user{id(i)},T,genDynVar,pararg{:});
    else
        out(i)=0;
    end
    if ~isempty(lastwarn)
        failed = [failed i];
    end
end
%---------------------------------------------------------

%------------------------------------------------------------
function [S,L] = singmat

% 0: testfunction must vanish
% 1: testfunction must not vanish
% everything else: ignore this testfunction

%   S = [  0 8 8
%          8 0 8
%          1 8 0 ];
%
%   L = [ 'BP'; 'H '; 'LP' ];
S = [ 0 8
    1 0];

L = [ 'BP'; 'LP' ];


function varargout = init(varargin)
x = varargin{1};
v = varargin{2};
WorkspaceInit(x,v);

% all done succesfully
varargout{1} = 0;
%---------------------------------------------------------
function done
WorkspaceDone;
%----------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ---------------------------------------------------------------

function [x,p] = rearr(x0)
%
% [x,p] = rearr(x0)
%
% Rearranges x0 into coordinates (x) and parameters (p)
global cds eds
nap = length(eds.ActiveParams);
ncoo = cds.ndim-nap;

p = eds.P0;
p(eds.ActiveParams) = x0((ncoo+1):end);
x = x0(1:ncoo);

function [x,p,T,initstates] = rearrnew(x0)
%
% [x,p] = rearr(x0)
%
% Rearranges x0 into coordinates (x) and parameters (p)
global cds eds ics
nap = length(eds.ActiveParams);
ncoo = cds.ndim-nap;

p = eds.P0;
p(eds.ActiveParams) = x0((ncoo+1):end);
initstates=p(ics.ParameterInitialStateCoordinates);
if ics.TimeHorizonContinuation
    p([ics.ParameterInitialStateCoordinates ics.ParameterEndTimeCoordinate])=[];
    T=x0((ncoo+1):end);
else
    T=p(end);
    p([ics.ParameterInitialStateCoordinates ics.ParameterEndTimeCoordinate])=[];
end
x = x0(1:ncoo);
% ---------------------------------------------------------

function WorkspaceInit(x,v)
global cds eds

% calculate some matrices to efficiently compute bialternate products (without loops)
% n = cds.ndim-1;
% a = reshape(1:(n^2),n,n);
% [bia,bin,bip] = bialt(a);
% if any(any(bip))
%     [eds.BiAlt_M1_I,eds.BiAlt_M1_J,eds.BiAlt_M1_V] = find(bip);
% else
%     eds.BiAlt_M1_I=1;eds.BiAlt_M1_J=1;eds.BiAlt_M1_V=n^2+1;
% end
% if any(any(bin))
%     [eds.BiAlt_M2_I,eds.BiAlt_M2_J,eds.BiAlt_M2_V] = find(bin);
% else
%      eds.BiAlt_M2_I=1;eds.BiAlt_M2_J=1;eds.BiAlt_M2_V=n^2+1;
% end
% if any(any(bia))
%     [eds.BiAlt_M3_I,eds.BiAlt_M3_J,eds.BiAlt_M3_V] = find(bia);
% else
%     eds.BiAlt_M3_I=1;eds.BiAlt_M3_J=1;eds.BiAlt_M3_V=n^2+1;
% end

% ------------------------------------------------------

function WorkspaceDone

% -------------------------------------------------------

function jac = DS_jac(DS_func,x)
global cds ics GlobalImpulseVariable
arcarg=GlobalImpulseVariable.arcarg;
arcnum=GlobalImpulseVariable.arcnum;
jumparg=GlobalImpulseVariable.jumparg;

jac=ics.sparseJac;
%%%%%
[x,p,T,initialstates] = rearrnew(arg);
%p = num2cell(p);
dynVar=x(ics.genDepVarCoordinates);
timepoints=[GlobalImpulseVariable.InitialTime x(GlobalImpulseVariable.JumpTimeCoordinates).' T];
for ii=1:arcnum+1
    dynVarLR=dynVar(:,(2*ii-1):2*ii); % variable representing the left and right limits
    if ii==1
        %Omega_0
        %         connecres=ics.jump(timepoints(ii),dynVarLR,p,jumparg(ii));
        %         initres=ics.initialcond(timepoints(ii),dynVarLR,p,arcarg(1),initialstates);
        res=ics.jump(timepoints(ii),dynVarLR,p,jumparg(ii));
        res=[res; ...
            ics.initialcond(timepoints(ii),dynVarLR,p,arcarg(1),initialstates)];
    elseif ii==arcnum+1
        %Omega_N
        %         connecres=[connecres;ics.jump(timepoints(ii),dynVarLR,p,jumparg(ii))];
        %         transres=ics.transversality(timepoints(ii),dynVarLR,p,arcarg(end));
        res=[res; ...
            ics.jump(timepoints(ii),dynVarLR,p,jumparg(ii))];
        res=[res;ics.transversality(timepoints(ii),dynVarLR,p,arcarg(end))];
    end
    if arcnum>1 && ii>=2 && ii<=arcnum
        %Omega_i
        %         connecres=[connecres;ics.jump(timepoints(ii),dynVarLR,p,jumparg(ii))];
        %         interiorjumpres=[interiorjumpres;ics.interiorjumpcond(timepoints(ii),dynVarLR,p,arcarg(ii-1),jumparg(ii))];
        res=[res; ...
            ics.jump(timepoints(ii),dynVarLR,p,jumparg(ii))];
        res=[res; ...
            ics.interiorjumpcond(timepoints(ii),dynVarLR,p,arcarg(ii-1),jumparg(ii))];
    end
    if ii<=arcnum
        %Gamma_i
        dynVarI=dynVar(:,2*ii:(2*ii+1)); % variable representing the points on the interval between two jumps
        %         dynres=[dynres;ics.continousdyn(timepoints(ii:ii+1),dynVarI,p,arcarg(ii))];
        res=[res;ics.continousdyn(timepoints(ii:ii+1),dynVarI,p,arcarg(ii))];
    end
end

%--------------------------------------------------------
function [x,v] = locate(id, x1, v1, x2, v2)
switch id
    case 1
        [x,v] = locateBP(id, x1, v1, x2, v2);
    otherwise
        error('No locator defined for singularity %d', id);
end

%----------------------------------------------------------------
function [out, failed] = testf(id, x, v)
global cds eds
ndim = cds.ndim;

if any(ismember(id,[1 2]))
    J = cjac(cds.curve_func,cds.curve_jacobian,x,[]);
end

out(3) = 0;
failed = [];

for i=id
    lastwarn('');

    switch i
        case 1 % BP
            % Jacobian extended with bordering vectors v and w
            B = [J; v'];
            out(1) = det(B);
            %
            %   case 2 % H
            %     A=J(:,1:ndim-1);
            %     A(:,end+1)=0;
            %     A1=sparse(eds.BiAlt_M1_I,eds.BiAlt_M1_J,A(eds.BiAlt_M1_V));
            %     A2=sparse(eds.BiAlt_M2_I,eds.BiAlt_M2_J,A(eds.BiAlt_M2_V));
            %     A3=sparse(eds.BiAlt_M3_I,eds.BiAlt_M3_J,A(eds.BiAlt_M3_V));
            %     A = A1-A2+A3;
            %     bigmat = [A eds.bigW; eds.bigV' eds.bigD];
            %
            %     Xg = bigmat \ [zeros(eds.nphase*(eds.nphase-1)/2,1); 1];
            %     out(2) = Xg(end);

        case 2%3 % LP
            out(2) = v(end);

        otherwise
            error('No such testfunction');
    end
    if ~isempty(lastwarn)
        failed = [failed i];
    end

end

%---------------------------------------------------------------------
function [failed,s] = process(id, x, v, s)
global cds
ndim = cds.ndim;

% WM: Removed SL array
fprintf('label = %s, x = ', s.label); printv(x);

switch id
    case 1 % BP
        s.data.v = v;
        s.msg  = sprintf('Branch point');
        %   case 2 % H
        %     s.data.lyapunov = lyapunov(x);
        %     if strcmp(s.data.lyapunov,'Neutral saddle')
        %         s.msg  = sprintf('Neutral saddle');
        %     else
        %         s.msg  = sprintf('Hopf');
        %         fprintf('First Lyapunov coefficient = %d\n', s.data.lyapunov);
        %     end
    case 2 % LP
        s.data.a = [];%a_lp(x);
        fprintf('a=%d\n',s.data.a);
        s.msg  = sprintf('Limit point');
end
% Compute eigenvalues for every singularity
J=cjac(cds.curve_func,cds.curve_jacobian,x,[]);
if ~issparse(J)
    [v,d]=eig(J(:,1:ndim-1));
else
    opt.disp=0;
    % WM: fixed a bug (incompatability between MatLab 6.0 and 5.5?)
    [v,d]=eigs(J(:,1:ndim-1),min(6,ndim-1),'lm',opt);
end

s.data.evec = v;
s.data.eval = diag(d)';

failed = 0;

% ---------------------------------------------------------------
function [x,v] = locateBP(id, x1, v1, x2, v2)
global cds

ndim = cds.ndim;
J = cjac(cds.curve_func,cds.curve_jacobian,x1,[]);
if ~issparse(J)
    [v,d]=eig(J(:,1:ndim-1)');
else
    opt.disp=0;
    [v,d]=eigs(J(:,1:ndim-1)', 'SM', opt);
end
[y,i]=min(abs(diag(d)));
p = v(:,i);
b = 0;
x = 0.5*(x1+x2);
i = 0;

u = [x; b; p];

[A,f]=locjac(x,b,p);
while i < cds.options.MaxCorrIters

    du = A\f;
    u = u - du;

    x = u(1:ndim);
    b = u(ndim+1);
    p = u(ndim+2:2*ndim);

    [A,f]=locjac(x,b,p);

    % WM: VarTol and FunTol were switched
    if norm(du) < cds.options.VarTolerance & norm(f) < cds.options.FunTolerance
        v = 0.5*(v1+v2);
        return;
    end
    i = i+1;
end
x=[];


%-----------------------------------------------------------------
function failed=saveintermediate(sout,xout,vout,hout,fout,contnum)
global eds cds ics
failed=0;
MODELINFO.OCMATICCONT=ics;
try
    if contnum==1
        save([ics.basicglobalvarfilename '4finiteimpulse'],'MODELINFO')
    end
    save([ics.basicresultfilename '4finiteimpulse'],'sout','xout','vout')
catch
    failed=1;
end
