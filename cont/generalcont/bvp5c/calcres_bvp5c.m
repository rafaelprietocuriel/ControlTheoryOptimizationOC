function RHS=calcres_bvp5c(x,y,freepar,modelpar,ode,bc)
%CALCRES_BVP4C  Evaluate the system of collocation equations res(y).
%   The derivative approximated at the midpoints and returned in Fmid is
%   used to estimate the residual.

% This function bases on bvp4c by
%   Jacek Kierzenka and Lawrence F. Shampine
%   Copyright 1984-2007 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $  $Date: 2007/05/23 18:54:08 $

global OCMATCONT OCBVP
nODEeval=0;
ya = y(:,1);
yb = y(:,end);
RHSbc = bc(ya,yb,freepar,modelpar);
N=OCBVP.N;
h=diff(x(1:OCBVP.nstages:end));
RHSode=zeros(OCBVP.neqn,N-1);

if OCMATCONT.OPTIONS.xyvectorized
    F = ode(x,y,OCMATCONT.HE.arcindex,freepar,modelpar);
    nODEeval = nODEeval + 1; % stats
else
    F = zeros(size(y));
    for i = 1 : N
        F(:,i) = ode(x(i),y(:,i),OCMATCONT.HE.arcindex,freepar,modelpar);
    end
    nODEeval = nODEeval + N; % stats
end
idx = 0;
for i = 1 : numel(h)
    ynn = y(:, idx + ones(1,OCBVP.nstages));
    ync = y(:, idx+2 : idx+OCBVP.nstages+1);

    hi = h(i);
    hA = hi*OCBVP.A;
    K = F(:, idx+1 : idx+OCBVP.nstages+1);

    % Form the residual in the Runge-Kutta formulas (collocation
    % equations) in terms of the intermediate solution values.
    RHSode(:, idx+1 : idx+OCBVP.nstages) = ynn + K*hA' - ync;

    idx = idx + OCBVP.nstages;
end

RHS = [RHSbc(:); RHSode(:)];

% colloc_RHS
OCBVP.F=F;