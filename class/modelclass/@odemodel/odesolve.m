    function [ocTrj,sol]=odesolve(odeObj,tspan,y0,opt,varargin)
%
% overloaded function to evaluate function which start with the modelname
% of odeObj.

arcid=[];
if isempty(odeObj)
    return
end
if nargin<4 || isempty(opt) 
    opt=defaultocoptions;
end
if nargin>=5
    arcid=varargin{1};
end
if isempty(arcid)
    arcid=0;
end
ODESolver=str2func(getocoptions(opt,'GENERAL','ODESolver'));
ODEOptions=opt.ODE;
if strcmp(ODEOptions.Events,'on')
    ODEOptions.Events=str2func([modelname(odeObj) 'Events']);
else
    ODEOptions.Events=[];
end
dynfunch=str2func([modelname(odeObj) 'Dynamics']);
par=parametervalue(odeObj);
sol=ODESolver(dynfunch,tspan,y0,ODEOptions,par,arcid);

sol.solverinfo.ODESolveroutput=sol;
sol.arcarg=arcid;
sol.arcinterval=sol.x([1 end]);
sol.arcposition=[1;length(sol.x)];
sol.x=sol.x/sol.x(end);
sol.modelparameter=par;
sol.modelname=modelname(odeObj);
ocTrj=octrajectory(sol);
