function [ocTrj,sol]=odesolve(ocObj,tspan,y0,opt,arcid,varargin)
%

dynfunc=[];
eventfunc=[];
if isempty(ocObj)
    ocTrj=octrajectory();
    return
end
if nargin<=3
    opt=defaultocoptions;
end
if nargin<=4
    arcid=0;
end
if nargin>=6
    dynfunc=varargin{1};
end
if nargin>=7
    eventfunc=varargin{2};
end
if isempty(dynfunc)
    dynfunch=str2func([modelname(ocObj) 'CanonicalSystem']);
end
if isempty(opt)
    opt=defaultocoptions;
end
if ischar(dynfunc)
    dynfunch=str2func([modelname(ocObj) dynfunc]);
end
ODESolver=str2func(getocoptions(opt,'GENERAL','ODESolver'));

ODEOptions=opt.ODE;
if strcmp(ODEOptions.Events,'on')
    if isempty(eventfunc)
        ODEOptions.Events=str2func([modelname(ocObj) 'Events']);
    else
        ODEOptions.Events=str2func([modelname(ocObj) eventfunc]);
    end
end
par=parametervalue(ocObj);
sol=ODESolver(dynfunch,tspan,y0,ODEOptions,par,arcid);

ocTrj=sol;
ocTrj.arcarg=arcid;
ocTrj.arcinterval=ocTrj.x([1 end]);
if ocTrj.arcinterval(1)>ocTrj.arcinterval(2)
    ocTrj.arcinterval=sort(ocTrj.arcinterval);
    ocTrj.x=ocTrj.x(end:-1:1);
    ocTrj.y=ocTrj.y(:,end:-1:1);
    if strcmp(getocoptions(opt,'GENERAL','ODESolver'),'ode45')
        sol.idata.f3d=sol.idata.f3d(:,:,end:-1:1);
    end
end
ocTrj.x=(ocTrj.x-ocTrj.x(1))/(ocTrj.x(end)-ocTrj.x(1));
ocTrj.arcposition=[1;length(ocTrj.x)];
ocTrj.userinfo=func2str(dynfunch);
ocTrj=octrajectory(ocTrj);
