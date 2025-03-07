function [dFdy,nFcalls,dFdp]=Fnumjac(ode,odeArgs)
%FNUMJAC  Numerically compute dF/dy and dF/dpar.
global OCBVP
% [dFdy,OCBVP.Joptions.fac,ignored,nFcalls] = ...
%     odenumjac(ode,odeArgs,OCBVP.Fref,OCBVP.Joptions);
nFcalls=[];
%OCBVP.Joptions.vectvars=[];
[dFdy,dum,dFdy_nfcn] = numjaccsd(ode,odeArgs,OCBVP.neqn,OCBVP.Joptions);
if OCBVP.explicitparameterdependence
%     [dFdp,OCBVP.dPoptions.fac] = ...
%         odenumjac(ode,odeArgs,OCBVP.Fref,OCBVP.dPoptions);
        [dFdp,dum,dFdp_nfcn] = numjaccsd(ode,odeArgs,OCBVP.neqn,OCBVP.dPoptions);
else
    dFdp=[];
    dFdp_nfcn=0;
end
% nFcalls = dFdy_nfcn + dFdp_nfcn;
