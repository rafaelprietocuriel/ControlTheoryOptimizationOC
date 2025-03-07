function out = ppdemodelequilibrium
%
% Standard ode file for MATCONT
                                                                                                  
out{1} = @init;
out{2} = @fun_eval;
out{3} = [];%@jacobian;
out{4} = [];%@jacobianp;
out{5} = [];%@hess;
out{6} = [];%@hessp;
out{7} = [];%@der3;
out{8} = [];%@der4;
out{9} = [];%@der5;
out{10}= @testadmissibility;
out{11}= @userfunc;
out{20}= @plotcont;
%
% --------------------------------------------------------------------------
function dxdt = fun_eval(t,depvar,varargin)
global OCMATFINITCONT
dxdt=OCMATFINITCONT.dynamics(depvar,[varargin{:}],OCMATFINITCONT.arcarg,OCMATFINITCONT.femdata);
% if any(abs(imag(dxdt)))
%     dxdt=repmat(nan,length(dxdt),1);
% end
% --------------------------------------------------------------------------
function [tspan,y0,options] = init
tspan = [0 10];
                                                                                                  
% --------------------------------------------------------------------------
function jac = jacobian(t,depvar,varargin)
global OCMATFINITCONT
jac=OCMATFINITCONT.jacobian(depvar,[varargin{:}],OCMATFINITCONT.arcarg,OCMATFINITCONT.femdata);
% --------------------------------------------------------------------------
function jacp = jacobianp(t,depvar,varargin)
global OCMATFINITCONT
jacp=OCMATFINITCONT.parameterjacobian(depvar,[varargin{:}],OCMATFINITCONT.arcarg);

% --------------------------------------------------------------------------
function h = hess(t,depvar,varargin)
global OCMATFINITCONT
h=OCMATFINITCONT.hessian(t,depvar,[varargin{:}],OCMATFINITCONT.arcarg);

                                                                                                  
% --------------------------------------------------------------------------
function hp = hessp(t,depvar,varargin)
global OCMATFINITCONT
hp=OCMATFINITCONT.parameterhessian(depvar,[varargin{:}],OCMATFINITCONT.arcarg);


% --------------------------------------------------------------------------
function h = der3(t,depvar,varargin)
global OCMATFINITCONT
h=OCMATFINITCONT.der3(t,depvar,[varargin{:}],OCMATFINITCONT.arcarg);

% --------------------------------------------------------------------------
function h = der4(t,depvar,varargin)
global OCMATFINITCONT
h=OCMATFINITCONT.der4(t,depvar,[varargin{:}],OCMATFINITCONT.arcarg);

% --------------------------------------------------------------------------
function h = der5(t,depvar,varargin)
global OCMATFINITCONT
h=OCMATFINITCONT.der5(t,depvar,[varargin{:}],OCMATFINITCONT.arcarg);

function h=plotcont(depvar,v,contnum)
global cds eds OCMATFINITCONT

switch func2str(cds.curve)
    case 'equilibrium'
        try
            h=OCMATFINITCONT.plotcont(depvar(1:end-1,1:contnum),OCMATFINITCONT.parametervalue,OCMATFINITCONT.arcarg,depvar(end,1:contnum),v);
        catch
            h=plot(depvar(cds.ndim,1:contnum),depvar(OCMATFINITCONT.plotcoord,1:contnum));
        end
    case 'hopf'
        h=plot3(depvar(cds.ndim-2,1:contnum),depvar(cds.ndim-1,1:contnum),depvar(cds.ndim,1:contnum));
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
