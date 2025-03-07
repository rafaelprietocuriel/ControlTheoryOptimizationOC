function out = modellimitpoint
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
out{10}= [];
out{10}= @testadmissibility;
out{11}= @userfunc;
out{20}= @plotcont;
%
% --------------------------------------------------------------------------
function dxdt = fun_eval(t,depvar,varargin)
global OCMATFINITCONT
%dxdt=OCMATFINITCONT.dynamics(t,depvar,[varargin{:}],OCMATFINITCONT.arcarg);
dxdt=OCMATFINITCONT.dynamics(depvar,[varargin{:}],OCMATFINITCONT.arcarg);
% --------------------------------------------------------------------------
function [tspan,y0,options] = init
tspan = [0 10];
                                                                                                  
% --------------------------------------------------------------------------
function jac = jacobian(t,depvar,varargin)
global OCMATFINITCONT
jac=OCMATFINITCONT.jacobian(depvar,[varargin{:}],OCMATFINITCONT.arcarg);

                                                                                                  
% --------------------------------------------------------------------------
function jacp = jacobianp(t,depvar,varargin)
global OCMATFINITCONT
jacp=OCMATFINITCONT.parameterjacobian(depvar,[varargin{:}],OCMATFINITCONT.arcarg);

function h=plotcont(x,v,contnum)
global cds lpds

h=plot3(x(cds.ndim-1,1:contnum),x(cds.ndim,1:contnum),x(1,1:contnum));

%----------------------------------------------------------------
function val=testadmissibility(t,depvar,varargin)
global OCMATFINITCONT
val=0;

try
    val=min(OCMATFINITCONT.testadmissibility(t,depvar,[varargin{:}],OCMATFINITCONT.arcarg));
end

%----------------------------------------------------------------
function val=userfunc(t,depvar,varargin)
global USERGLOBAL lpds
val=0;

try
    if isfield(USERGLOBAL,'UserValue')
        if isfield(USERGLOBAL,'UserValueIndex')
            val=(USERGLOBAL.UserValue-varargin{lpds.ActiveParams(USERGLOBAL.UserValueIndex)});
        else
            val=(USERGLOBAL.UserValue-varargin{lpds.ActiveParams(1)})*(USERGLOBAL.UserValue-varargin{lpds.ActiveParams(2)});
        end
    elseif isfield(USERGLOBAL,'UserValue1') && isfield(USERGLOBAL,'UserValue2')
        val=(USERGLOBAL.UserValue1-varargin{lpds.ActiveParams(1)})*(USERGLOBAL.UserValue2-varargin{lpds.ActiveParams(1)})*(USERGLOBAL.UserValue1-varargin{lpds.ActiveParams(2)})*(USERGLOBAL.UserValue2-varargin{lpds.ActiveParams(2)});
    end

end
