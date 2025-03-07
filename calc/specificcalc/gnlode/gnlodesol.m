function out=gnlodesol()
%
	
out{1}=@statedynamics;
out{2}=@adjointdynamics;
out{3}=@gradienthamiltonian;
out{4}=@objectivefunction;
out{5}=@initialstate;
out{6}=@transversalitycondition;
out{7}=@integralconstraint;
out{8}=@adjointintegral;
	
out{11}=@plotcalcoc;
out{13}=@saveintermediate;
out{14}=@drearr;
out{22}=@formatsolution;
	
%--------------------------------------------------------------------------
function [dnlxdt dlxdt]=statedynamics(x,t,y,z,p,q,u,v,w,varargin)
global OCMATGNLODE

[nldepvar,ldepvar,nlintdepvar,lintdepvar,nlctrl,lctrl,w,modelpar,exp_rt]= ...
    drearr(y,z,p,q,[],[],[],[],u,v,w);

[dnlxdt dlxdt]=OCMATGNLODE.statedynamics(x,t,nldepvar,ldepvar,nlintdepvar,lintdepvar,nlctrl,lctrl,modelpar);
	
%--------------------------------------------------------------------------
function [dnlldt dlldt]=adjointdynamics(x,t,y,z,p,q,xi,theta,eta,zeta,u,v,w,varargin)
global OCMATGNLODE
[nldepvar,ldepvar,nlintdepvar,lintdepvar,nlctrl,lctrl,w,modelpar,exp_rt]= ...
    drearr(y,z,p,q,xi,theta,eta,zeta,u,v,w);

[dnlldt dlldt]=OCMATGNLODE.adjointdynamics(x,t,nldepvar,ldepvar,nlintdepvar,lintdepvar,nlctrl,lctrl,modelpar);
	
%--------------------------------------------------------------------------
function [eta zeta]=adjointintegral(x,t,y,z,p,q,xi,theta,u,v,w,varargin)
global OCMATGNLODE
[nldepvar,ldepvar,nlintdepvar,lintdepvar,nlctrl,lctrl,w,modelpar]= ...
    drearr(y,z,p,q,xi,theta,[],[],u,v,w);

[eta zeta]=OCMATGNLODE.adjointintegral(x,t,nldepvar,ldepvar,nlintdepvar,lintdepvar,nlctrl,lctrl,modelpar);
	
%--------------------------------------------------------------------------
function [dHdu,dHdv,dHdw]=gradienthamiltonian(x,t,y,z,p,q,xi,theta,eta,zeta,u,v,w,varargin)
global OCMATGNLODE
[nldepvar,ldepvar,nlintdepvar,lintdepvar,nlctrl,lctrl,w,modelpar,exp_rt]= ...
    drearr(y,z,p,q,xi,theta,eta,zeta,u,v,w);
[dHdu,dHdv,dHdw]=OCMATGNLODE.gradienthamiltonian(x,t,nldepvar,ldepvar,nlintdepvar,lintdepvar,nlctrl,lctrl,modelpar);

%--------------------------------------------------------------------------
function J=objectivefunction(x,t,y,z,p,q,u,v,w,varargin)
global OCMATGNLODE
[nldepvar,ldepvar,nlintdepvar,lintdepvar,nlctrl,lctrl,w,modelpar,exp_rt]= ...
    drearr(y,z,p,q,[],[],[],[],u,v,w);

B=OCMATGNLODE.objectivefunction(x,t,nldepvar,ldepvar,nlintdepvar,lintdepvar,nlctrl,lctrl,modelpar).';
J=trapzint1q(OCMATGNLODE.dt(:),exp_rt(:).*B(:));

%--------------------------------------------------------------------------
function z0=initialstate(w,varargin)
global OCMATGNLODE
z0=[];
if strcmp(OCMATGNLODE.initialstate,'fix')
    z0=OCMATGNLODE.z0;
end

%--------------------------------------------------------------------------
function tc=transversalitycondition(x,y,z,p,q,u,v,w,varargin)
global OCMATGNLODE OCMATCALC
[nldepvar,ldepvar,nlintdepvar,lintdepvar,nlctrl,lctrl,w,modelpar,exp_rt]= ...
    drearr(y,z,p,q,[],[],[],[],u,v,w);

tc=OCMATGNLODE.transversalitycondition(OCMATCALC.tspan,depvar,ctrl,modelpar);


%--------------------------------------------------------------------------
function [nlint lint]=integralconstraint(x,t,y,z,u,v,w,varargin)
global OCMATGNLODE

[nldepvar,ldepvar,nlintdepvar,lintdepvar,nlctrl,lctrl,w,modelpar,exp_rt]= ...
    drearr(y,z,[],[],[],[],[],[],u,v,w);

[nlint lint]=OCMATGNLODE.integralconstraint(x,t,nldepvar,ldepvar,nlintdepvar,lintdepvar,nlctrl,lctrl,modelpar);
	
%--------------------------------------------------------------------------
function h=plotcalcoc(x,t,y,z,p,q,xi,theta,eta,zeta,u,v,w,Val,du,dv,dw)
global OCMATGNLODE
[nldepvar,ldepvar,nlintdepvar,lintdepvar,nlctrl,lctrl,w,modelpar]= ...
    drearr(y,z,p,q,xi,theta,eta,zeta,u,v,w);

h=OCMATGNLODE.plotcalcoc(x,t,nldepvar,ldepvar,nlintdepvar,lintdepvar,nlctrl,lctrl,modelpar,Val,du,dv,dw);
figure(gcf)

%-----------------------------------------------------------------
function failed=saveintermediate(sout)
global OCMATGNLODE
failed=0;
try
    save([OCMATGNLODE.basicresultfilename '4godesol'],'sout')
catch
    failed=1;
end


% ------------------------------------------------------
function [nldepvar,ldepvar,nlintdepvar,lintdepvar,nlctrl,lctrl,w,modelpar,exp_rt]=drearr(y,z,p,q,xi,theta,eta,zeta,u,v,w)
global OCMATGNLODE 
modelpar=OCMATGNLODE.modelparameter;
exp_rt=OCMATGNLODE.exp_rt;
nlctrl=u;
lctrl=v;
if isempty(y)
    nldepvar=[];
else
    nldepvar=[y xi];
end
if isempty(z)
    ldepvar=[];
else
    ldepvar=[z theta];
end
if isempty(p)
    nlintdepvar=[];
else
    nlintdepvar=[p eta];
end
if isempty(q)
    lintdepvar=[];
else
    lintdepvar=[q zeta];
end




%----------------------------------------------------------------
function out=formatsolution(x,t,y,z,p,q,xi,theta,eta,zeta,u,v,w,Val,du,dv,dw)
global OCMATGNLODE OCMATCALC
[nldepvar,ldepvar,nlintdepvar,lintdepvar,nlctrl,lctrl,w,modelpar]= ...
    drearr(y,z,p,q,xi,theta,eta,zeta,u,v,w);

out.t=t;
out.x=x;
out.y=nldepvar;
out.z=ldepvar;
out.p=nlintdepvar;
out.q=lintdepvar;
out.u=nlctrl;
out.v=lctrl;
out.w=w;
out.modelinfo.modelname=OCMATCALC.modelname;
out.modelinfo.modelparameter=modelpar;
out.modelinfo.modeltype=OCMATCALC.modeltype;
out.modelinfo.optimizationtype=OCMATCALC.optimizationtype;
if ~isempty(out.y)
    out.modelinfo.y0=out.y(:,1:end/2,1);
end
if ~isempty(out.z)
    out.modelinfo.z0=out.z(1:end/2,1);
end

out.numericalinfo.du=du;
out.numericalinfo.dv=dv;
out.numericalinfo.dw=dw;
out.numericalinfo.objectivevalue=Val;

out.solverinfo.gradmap=func2str(OCMATCALC.gradmap);
out.solverinfo.problemfunction=OCMATCALC.problem_func;
out.solverinfo.opt=OCMATCALC;
out.solverinfo.solver=func2str(OCMATCALC.GradSolver);

