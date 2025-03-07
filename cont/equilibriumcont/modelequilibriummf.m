function out = modelequilibriummf
%
% Standard ode file for MATCONT
                                                                                                  
out{1} = @init;
out{2} = @fun_eval;
out{3} = @jacobian;
out{4} = @jacobianp;
out{5} = [];
out{6} = [];
out{7} = [];
out{8} = [];
out{9} = [];
out{10}= @testadmissibility;
out{20}= @plotcont;
%
% --------------------------------------------------------------------------
function dxdt = fun_eval(t,depvar,varargin)
global OCMATFINITCONT
[depvar par]=rearr(depvar,[varargin{:}]);
dxdt=OCMATFINITCONT.dynamics(depvar,par,OCMATFINITCONT.arcarg);
dxdt(OCMATFINITCONT.equilibriummffreecoord)=[];
% --------------------------------------------------------------------------
function [tspan,y0,options] = init
tspan = [0 10];
                                                                                                  
% --------------------------------------------------------------------------
function jac = jacobian(t,depvar,varargin)
global OCMATFINITCONT
[depvar par]=rearr(depvar,[varargin{:}]);
jac=OCMATFINITCONT.jacobian(depvar,par,OCMATFINITCONT.arcarg);
jac(OCMATFINITCONT.equilibriummffreecoord,:)=[];
jac(:,OCMATFINITCONT.equilibriummffreecoord)=[];
                                                                                                  
% --------------------------------------------------------------------------
function jacp = jacobianp(t,depvar,varargin)
global OCMATFINITCONT
[depvar par]=rearr(depvar,[varargin{:}]);
jac=OCMATFINITCONT.jacobian(depvar,par,OCMATFINITCONT.arcarg);
jacp=OCMATFINITCONT.parameterjacobian(depvar,par,OCMATFINITCONT.arcarg);
jacp=[jacp jac(:,OCMATFINITCONT.equilibriummffreecoord)];
jacp(OCMATFINITCONT.equilibriummffreecoord,:)=[];

function h=plotcont(x,v,contnum)
global cds eds

switch func2str(cds.curve)
    case 'equilibrium'
        h=plot(x(cds.ndim,1:contnum),x(1,1:contnum));
    case 'hopf'
        h=plot(x(cds.ndim-2,1:contnum),x(cds.ndim-1,1:contnum));
end

%----------------------------------------------------------------
function val=testadmissibility(t,depvar,varargin)
global OCMATFINITCONT
[depvar par]=rearr(depvar,[varargin{:}]);

val=min(OCMATFINITCONT.testadmissibility(t,depvar,par,OCMATFINITCONT.arcarg));

%----------------------------------------------------------------

function [depvar par]=rearr(depvar,totalpar)
global OCMATFINITCONT
depvar(OCMATFINITCONT.equilibriummfcoord)=depvar;
depvar(OCMATFINITCONT.equilibriummffreecoord)=totalpar(OCMATFINITCONT.equilibriummfcoordposition);
par=totalpar(OCMATFINITCONT.parametercoord);


