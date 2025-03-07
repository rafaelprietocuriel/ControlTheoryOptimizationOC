function out = odemodelequilibrium
%
% Standard ode file for MATCONT
                                                                                                  
out{1} = @init;
out{2} = @fun_eval;
out{3} = [];%@jacobian;
out{4} = [];%@jacobianp;
out{5} = [];
out{6} = [];
out{7} = [];
out{8} = [];
out{9} = [];
out{10}= @testadmissibility;
out{11}= @userfunc;
out{20}= @plotcont;
%
% --------------------------------------------------------------------------
function dxdt = fun_eval(t,depvar,varargin)
global OCMATFINITCONT
dxdt=OCMATFINITCONT.dynamics(t,depvar,[varargin{:}],OCMATFINITCONT.arcarg);
% --------------------------------------------------------------------------
function [tspan,y0,options] = init
tspan = [0 10];
                                                                                                  
% --------------------------------------------------------------------------
function jac = jacobian(t,depvar,varargin)
global OCMATFINITCONT
jac=OCMATFINITCONT.jacobian(t,depvar,[varargin{:}],OCMATFINITCONT.arcarg);

                                                                                                  
% --------------------------------------------------------------------------
function jacp = jacobianp(t,depvar,varargin)
global OCMATFINITCONT
jacp=OCMATFINITCONT.parameterjacobian(t,depvar,[varargin{:}],OCMATFINITCONT.arcarg);

function h=plotcont(depvar,v,contnum)
global cds eds OCMATFINITCONT

switch func2str(cds.curve)
    case 'equilibrium'
        try
            h=OCMATFINITCONT.plotcont(depvar(1:cds.ndim-1,1:contnum),OCMATFINITCONT.parametervalue,OCMATFINITCONT.arcarg,depvar(end,1:contnum),v);
        catch
            h=plot(depvar(cds.ndim,1:contnum),depvar(OCMATFINITCONT.plotcoord,1:contnum));
        end
    case 'hopf'
        h=plot(depvar(cds.ndim-2,1:contnum),x(cds.ndim-1,1:contnum));
end

%----------------------------------------------------------------
function val=testadmissibility(t,depvar,varargin)
global OCMATFINITCONT
val=0;

try
    val=min(OCMATFINITCONT.testadmissibility(t,depvar,[varargin{:}],OCMATFINITCONT.arcarg));
end
%----------------------------------------------------------------
function val=userfunc(t,depvar,varargin)
global USERGLOBAL eds
val=0;

try
    if isfield(USERGLOBAL,'UserValue')
        val=USERGLOBAL.UserValue-varargin{eds.ActiveParams};
    elseif isfield(USERGLOBAL,'UserValue1') && isfield(USERGLOBAL,'UserValue2')
        val=(USERGLOBAL.UserValue1-varargin{eds.ActiveParams})*(USERGLOBAL.UserValue2-varargin{eds.ActiveParams});
    end

end