function out=godesol()
%
	
out{1}=@statedynamics;
out{2}=@adjointdynamics;
out{3}=@gradienthamiltonian;
out{4}=@objectivefunction;
out{5}=@initialstate;
out{6}=@transversalitycondition;
	
out{11}=@plotcalcoc;
out{13}=@saveintermediate;
out{14}=@drearr;
	
%--------------------------------------------------------------------------
function dxdt=statedynamics(t,z,v,varargin)
global OCMATGODE

[depvar,ctrl,w,modelpar]=drearr(z,[],v,[]);

dxdt=OCMATGODE.statedynamics(t,depvar,ctrl,modelpar);
	
%--------------------------------------------------------------------------
function dldt=adjointdynamics(t,z,theta,v,w,varargin)
global OCMATGODE
[depvar,ctrl,w,modelpar]=drearr(z,theta,v,w);

dldt=OCMATGODE.adjointdynamics(t,depvar,ctrl,modelpar);
	
%--------------------------------------------------------------------------
function [dHdu,dHdv,dHdw]=gradienthamiltonian(t,z,theta,v,w,varargin)
global OCMATGODE
[depvar,ctrl,w,modelpar]=drearr(z,theta,v,w);
dHdu=[];
dHdw=[];
dHdv=OCMATGODE.gradienthamiltonian(t,depvar,ctrl,modelpar);
	
%--------------------------------------------------------------------------
function J=objectivefunction(t,z,v,w,varargin)
global OCMATGODE
[depvar,ctrl,w,modelpar,exp_rt]=drearr(z,[],v,w);

B=OCMATGODE.objectivefunction(t,depvar,ctrl,modelpar).';
J=trapzint1q(OCMATGODE.dt(:),exp_rt(:).*B(:));	

%--------------------------------------------------------------------------
function z0=initialstate(w,varargin)
global OCMATGODE
z0=[];
if strcmp(OCMATGODE.initialstate,'fix')
    z0=OCMATGODE.z0;
end

%--------------------------------------------------------------------------
function tc=transversalitycondition(z,v,w,varargin)
global OCMATGODE OCMATCALC
[depvar,ctrl,w,modelpar]=drearr(z,[],v,w);

tc=OCMATGODE.transversalitycondition(OCMATCALC.tspan,depvar,ctrl,modelpar);
	
	
%--------------------------------------------------------------------------
function h=plotcalcoc(t,u,v,w,Val,du,dv,dw,z,theta)
global OCMATGODE
[depvar,ctrl,w,modelpar]=drearr(z,theta,v,w);

h=OCMATGODE.plotcalcoc(t,depvar,ctrl,modelpar,Val,du,dv,dw);
figure(gcf)

%-----------------------------------------------------------------
function failed=saveintermediate(sout)
global OCMATGODE
failed=0;
try
    save([OCMATGODE.basicresultfilename '4godesol'],'sout')
catch
    failed=1;
end


% ------------------------------------------------------
function [depvar,ctrl,w,modelpar,exp_rt]=drearr(z,theta,v,w)
global OCMATGODE 
modelpar=OCMATGODE.modelparameter;
exp_rt=OCMATGODE.exp_rt;
ctrl=v;
depvar=[z;theta];
