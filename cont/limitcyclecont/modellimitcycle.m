function out = modellimitcycle()
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
out{11}= @userfunc;
out{20}= @plotcont;
%
% --------------------------------------------------------------------------
function dxdt = fun_eval(t,depvar,varargin)
global OCMATFINITCONT
dxdt=OCMATFINITCONT.dynamics(t,depvar,[varargin{:}],OCMATFINITCONT.arcarg);
%dxdt=real(dxdt);
% --------------------------------------------------------------------------
function [tspan,y0,options] = init
tspan = [0 10];
                                                                                                  
% --------------------------------------------------------------------------
function jac = jacobian(t,depvar,varargin)
global OCMATFINITCONT
jac=OCMATFINITCONT.jacobian(depvar,[varargin{:}],OCMATFINITCONT.arcarg);
%jac=real(jac);

                                                                                                  
% --------------------------------------------------------------------------
function jacp = jacobianp(t,depvar,varargin)
global OCMATFINITCONT
jacp=OCMATFINITCONT.parameterjacobian(depvar,[varargin{:}],OCMATFINITCONT.arcarg);
%jacp=real(jacp);

function h=plotcont(x,v,contnum)
global cds lds

subplot(1,2,1)
h=plot(x(cds.ndim,1:contnum),x(cds.ndim-1,1:contnum));
subplot(1,2,2)
y=reshape(x(1:cds.ndim-2,contnum),lds.nphase,lds.tps);
h(2)=plot(y(1,:),y(2,:));
ylabel('T')
title(['Step : ' num2str(contnum)])

%----------------------------------------------------------------
function val=testadmissibility(t,depvar,varargin)
global OCMATFINITCONT
val=0;

try
    val=min(OCMATFINITCONT.testadmissibility(t,depvar,[varargin{:}],OCMATFINITCONT.arcarg));
end


%----------------------------------------------------------------
function val=userfunc(t,depvar,varargin)
global USERGLOBAL lds
val=0;

try
    if isfield(USERGLOBAL,'UserValue')
        val=USERGLOBAL.UserValue-varargin{lds.ActiveParams};
    elseif isfield(USERGLOBAL,'UserValue1') && isfield(USERGLOBAL,'UserValue2')
        val=(USERGLOBAL.UserValue1-varargin{lds.ActiveParams})*(USERGLOBAL.UserValue2-varargin{lds.ActiveParams});
    end

end
