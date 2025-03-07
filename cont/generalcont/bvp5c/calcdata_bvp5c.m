function [K1 ymid]=calcdata_bvp5c(X,Y,freepar,modelpar,ode)
%CALCRES_BVP4C  Evaluate the system of collocation equations res(y).
%   The derivative approximated at the midpoints and returned in Fmid is
%   used to estimate the residual.

global OCMATCONT OCBVP
N=numel(X);
nstages=OCBVP.nstages;
sq5=OCBVP.sq5;

if OCMATCONT.OPTIONS.xyvectorized
    F = ode(X,Y,OCMATCONT.HE.arcindex,freepar,modelpar);
else
    F = zeros(size(Y));
    for i = 1 : N
        F(:,i) = ode(X(i),Y(:,i),OCMATCONT.HE.arcindex,freepar,modelpar);
    end
end

h = diff(X(1:nstages:end));
y = Y(:,1:nstages:end);
K1 = F(:,1:nstages:end);
K2 = F(:,2:nstages:end);
K3 = F(:,3:nstages:end);
ymid = y(:,1:end-1) + ...
    ( 17/192*K1(:,1:end-1) + (40+15*sq5)/192*K2 + ...
    (40-15*sq5)/192*K3 - 1/192*K1(:,2:end)) * ...
    spdiags(h(:),0,numel(h),numel(h));
