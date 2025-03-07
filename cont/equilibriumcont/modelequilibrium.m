function out = modelequilibrium
%
% Standard ode file for MATCONT
                                                                                                  
out{1} = @init;
out{2} = @fun_eval;
out{3} = @jacobian;
out{4} = @jacobianp;
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
if isfield(OCMATFINITCONT,'combinationfunction') && ~isempty(OCMATFINITCONT.combinationfunction)
    par=[varargin{:}];
    par=OCMATFINITCONT.combinationfunction(par);
else
    par=[varargin{:}];
end
dxdt=OCMATFINITCONT.dynamics(depvar,par,OCMATFINITCONT.arcarg);
% if any(abs(imag(dxdt)))
%     dxdt=repmat(nan,length(dxdt),1);
% end
% --------------------------------------------------------------------------
function [tspan,y0,options] = init
tspan = [0 10];
                                                                                                  
% --------------------------------------------------------------------------
function jac = jacobian(t,depvar,varargin)
global OCMATFINITCONT
% numJacOpt.diffvar=1;
% numJacOpt.vectvars=[];
% jac=numjaccsd(OCMATFINITCONT.dynamics,{depvar,[varargin{:}],OCMATFINITCONT.arcarg},numel(depvar),numJacOpt);
if isfield(OCMATFINITCONT,'combinationfunction') && ~isempty(OCMATFINITCONT.combinationfunction)
    par=[varargin{:}];
    par=OCMATFINITCONT.combinationfunction(par);
else
    par=[varargin{:}];
end

jac=OCMATFINITCONT.jacobian(depvar,par,OCMATFINITCONT.arcarg);
% --------------------------------------------------------------------------
function jacp = jacobianp(t,depvar,varargin)
global OCMATFINITCONT
% numJacOpt.diffvar=2;
% numJacOpt.vectvars=[];
% jacp=numjaccsd(OCMATFINITCONT.dynamics,{depvar,[varargin{:}].',OCMATFINITCONT.arcarg},numel(depvar),numJacOpt);
if isfield(OCMATFINITCONT,'combinationfunction') && ~isempty(OCMATFINITCONT.combinationfunction)
    par=[varargin{:}];
    numJacOpt.diffvar=1;
    numJacOpt.vectvars=[];
    %dpfdmu=numjaccsd(OCMATFINITCONT.combinationfunction,{depvar,par.',OCMATFINITCONT.arcarg},numel(depvar),numJacOpt);
    dpfdmu=numjaccsd(OCMATFINITCONT.combinationfunction,{par.'},numel(par),numJacOpt);
else
    par=[varargin{:}];
    dpfdmu=1;
end

jacp=OCMATFINITCONT.parameterjacobian(depvar,par,OCMATFINITCONT.arcarg);
jacp=jacp*dpfdmu;
%jacp=jacp(:,OCMATFINITCONT.activeparameter);

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
        xlabel(OCMATFINITCONT.parametername{1},'Interpreter','latex')
    case 'hopf'
        %h=plot3(depvar(cds.ndim-1,1:contnum),depvar(cds.ndim-2,1:contnum),depvar(1,1:contnum));
        h=plot(depvar(cds.ndim-1,1:contnum),depvar(cds.ndim-2,1:contnum));
        xlabel(OCMATFINITCONT.parametername{2},'Interpreter','latex')
        ylabel(OCMATFINITCONT.parametername{1},'Interpreter','latex')
end

%----------------------------------------------------------------
function val=testadmissibility(t,depvar,varargin)
global OCMATFINITCONT
val=0;
if isfield(OCMATFINITCONT,'combinationfunction') && ~isempty(OCMATFINITCONT.combinationfunction)
    par=[varargin{:}];
    par=OCMATFINITCONT.combinationfunction(par);
else
    par=[varargin{:}];
end

try
    val=min(OCMATFINITCONT.testadmissibility(t,depvar,par,OCMATFINITCONT.arcarg))+1e-8;
end

%----------------------------------------------------------------
function val=userfunc(t,depvar,varargin)
global USERGLOBAL eds hds OCMATFINITCONT
%hds=[];
val=0;

try
    if isfield(OCMATFINITCONT,'conttype') && strcmp(OCMATFINITCONT.conttype,'hopf')
        if isfield(USERGLOBAL,'UserValue')
            if isfield(USERGLOBAL,'UserValueIndex')
                val=(USERGLOBAL.UserValue-varargin{hds.ActiveParams(USERGLOBAL.UserValueIndex)});
            else
                val=(USERGLOBAL.UserValue-varargin{hds.ActiveParams(1)})*(USERGLOBAL.UserValue-varargin{hds.ActiveParams(2)});
            end
        elseif isfield(USERGLOBAL,'UserValue1') && isfield(USERGLOBAL,'UserValue2')
            val=(USERGLOBAL.UserValue1-varargin{hds.ActiveParams(1)})*(USERGLOBAL.UserValue2-varargin{hds.ActiveParams(1)})*(USERGLOBAL.UserValue1-varargin{hds.ActiveParams(2)})*(USERGLOBAL.UserValue2-varargin{hds.ActiveParams(2)});
        end
    elseif isfield(OCMATFINITCONT,'conttype') && strcmp(OCMATFINITCONT.conttype,'equilibrium')
        if isfield(USERGLOBAL,'UserValue')
            val=USERGLOBAL.UserValue-varargin{eds.ActiveParams(1)};
        elseif isfield(USERGLOBAL,'UserValue1') && isfield(USERGLOBAL,'UserValue2')
            val=(USERGLOBAL.UserValue1-varargin{eds.ActiveParams})*(USERGLOBAL.UserValue2-varargin{eds.ActiveParams});
        end
    else
        if isfield(USERGLOBAL,'UserValue') && ~isempty(hds)
            val=USERGLOBAL.UserValue-varargin{hds.ActiveParams(1)};
        elseif isfield(USERGLOBAL,'UserValue') && ~isempty(eds)
            val=USERGLOBAL.UserValue-varargin{eds.ActiveParams};
        elseif isfield(USERGLOBAL,'UserValue1') && isfield(USERGLOBAL,'UserValue2') && ~isempty(hds)
            val=(USERGLOBAL.UserValue1-varargin{hds.ActiveParams(1)})*(USERGLOBAL.UserValue2-varargin{hds.ActiveParams(1)});
        elseif isfield(USERGLOBAL,'UserValue1') && isfield(USERGLOBAL,'UserValue2') && isempty(hds)
            val=(USERGLOBAL.UserValue1-varargin{eds.ActiveParams})*(USERGLOBAL.UserValue2-varargin{eds.ActiveParams});
        end
    end
end

