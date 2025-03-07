function prepare4meshadapt_bvp5c(X,Y)
global OCBVP

F=OCBVP.F;
sq5=OCBVP.sq5;
nstages=OCBVP.nstages;
h = diff(X(1:nstages:end));
y = Y(:,1:nstages:end);
K1 = F(:,1:nstages:end);
K2 = F(:,2:nstages:end);
K3 = F(:,3:nstages:end);
OCBVP.ymid = y(:,1:end-1) + ...
    ( 17/192*K1(:,1:end-1) + (40+15*sq5)/192*K2 + ...
    (40-15*sq5)/192*K3 - 1/192*K1(:,2:end)) * ...
    spdiags(h(:),0,numel(h),numel(h));
