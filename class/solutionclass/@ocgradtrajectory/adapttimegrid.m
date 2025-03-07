function [ocgTrjNew,tnew,errval]=adapttimegrid(ocgTrj,ocgTrjFine,varargin)
%
% [OCGTRJNEW,TNEW,ERRVAL]=ADAPTTIMEGRID(OCGTRJ,OCGTRJFINE)
%
% adapt the time mesh:
% ocgTrj is the base solution, ocgTrjFine is calculated on the refined time
% grid of ocgTrj (h -> h/2). The local error is evaluated
%       err_i:=||Y_i-Y^f_i||, 
% where Y consists of the state and costate values. Using the theory of
% adaptive step size for Runge Kutta schemes, we find at each time point
% the adapted time step
% ha_i=(tol/err_i)^(1/(p+1))*h_i, where p is the order of the Runge-Kutta
% scheme and tol the provided tolerance.
% The new time hn grid is calculated solving
%
% min (hn-ha_i)'*(hn-ha_i), s.t. sum(hn_i)=T, hn_i>=min_step
%
% As a second criterion it can be specified that the new time grid is
% equally distributed, which can be reformulated as
% min (hn-1/(n-1))'*(hn-1/(n-1)), s.t. sum(hn_i)=T, hn_i>=min_step
% 
% In total a weighted sum of both criteria can be formulated
% min (hn-rho*ha_i-(1-rho)/(n-1))'*(hn-rho*ha_i-(1-rho)/(n-1)), s.t.
% sum(hn_i)=T, hn_i>=min_step, 0<=rho<=1.

ocgTrjNew=[];
errval=[];
tnew=[];

tol=[];
mindist=[];
opt=[];
disp=[];
rho=[];
scale=[];
lambda=[];

if isempty(ocgTrj) || isempty(ocgTrjFine)
    return
end
odeslver=odesolver(ocgTrj);
t=time(ocgTrj);
h=diff(t);
tf=time(ocgTrjFine);
n=length(t);

if length(tf)~=2*n-1
    return
end
    
tolidx=find(strcmpi(varargin,'tol'));
scaleidx=find(strcmpi(varargin,'scale'));
lambdaidx=find(strcmpi(varargin,'lambda'));
dispidx=find(strcmpi(varargin,'disp'));
mindistidx=find(strcmpi(varargin,'mindist'));
optionidx=find(strcmpi(varargin,'option'));
rhoidx=find(strcmpi(varargin,'rho'));
if ~isempty(tolidx)
    tol=varargin{tolidx+1};
end
if ~isempty(rhoidx)
    rho=varargin{rhoidx+1};
end
if ~isempty(mindistidx)
    mindist=varargin{mindistidx+1};
end
if ~isempty(dispidx)
    disp=varargin{dispidx+1};
end
if ~isempty(lambdaidx)
    lambda=varargin{lambdaidx+1};
end
if ~isempty(scaleidx)
    scale=varargin{scaleidx+1};
end
if ~isempty(optionidx)
    opt=varargin{optionidx+1};
end
if isempty(tol)
    tol=1e-4;
end
if isempty(rho)
    rho=0.5;
end
if isempty(disp)
    disp=0;
end
if isempty(opt)
    %opt=defaultocoptions;
    opt=setocoptions('EQ','LargeScale','off');
end
if isempty(scale)
    scale='lin';
end
if strcmp(scale,'log')
    if isempty(lambda)
        return
    end
    lambda=-abs(real(lambda));
end
% X=[state(ocgTrj); ...
%     costate(ocgTrj); ...
%     control(ocgTrj)];
% 
% Xf=[state(ocgTrjFine); ...
%     costate(ocgTrjFine); ...
%     control(ocgTrjFine)];
X=[state(ocgTrj); ...
    costate(ocgTrj)];

Xf=[state(ocgTrjFine); ...
    costate(ocgTrjFine)];

switch scale
    case 'lin' 
        href=1/(n-1);
    case 'log'
        T=t(end)-t(1);
        XT=1e-5;
        X0=XT*exp(-lambda*T);
        href=diff([0 log(((n-1)-(1:(n-1)))/(n-1)+(1:(n-1))/(n-1)*XT/X0)/lambda]);
end
stepweigth=href;

errval=X(:,1:end)-Xf(:,1:2:end);

normerrval=sqrt(sum(errval.^2));
normerrval(1)=[];

switch odeslver
    case {'heunloc','ode2loc'}
        p=2;
    case 'ode5loc'
        p=5;
    otherwise
        p=1;
end
if isempty(mindist)
    mindist=(1e-15)^(1/(p+1));
end
errorweigth=(tol./normerrval).^(1/(p+1)).*h;

wheigth=-(rho*errorweigth+(1-rho)*stepweigth);
hnew=quadprog(1,wheigth.',[],[],ones(1,n-1),t(end),repmat(mindist,n-1,1),repmat(t(end),n-1,1),[],opt.EQ);
tnew=cumsum([t(1) hnew.']);
tnew(end)=t(end);
ocgTrjNew=deval(ocgTrj,tnew);

if disp
    plot(t(2:end),normerrval)
    %plot(tmid,normDX-meanlength)
    figure(gcf)
    pause
end