function res=residual_bvp5c(X,Y,freepar,modelpar,RHS,ode)
global OCBVP

%RESIDUAL_ESTIMATE  Estimate the residual in each subinterval

% Compute the max norm of the residual at Cres points in each subinterval.
% The norm is computed with the value of the solution (and thresh) as weights.
% When we trust the asymptotic behavior, the residual is only sampled at
% the midpoint. Otherwise, the residual is sampled at its three local max.
% ymid is the solution at midpoints.

% This function bases on bvp4c by
%   Jacek Kierzenka and Lawrence F. Shampine
%   Copyright 1984-2007 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $  $Date: 2007/05/23 18:54:08 $
solutionVerificationStage=OCBVP.solutionVerificationStage;
F=OCBVP.F;
ymid=OCBVP.ymid;
rtol=OCBVP.rtol;

if solutionVerificationStage
    % use three error samples per mesh interval
    nsamples = 3;
    res = residualEstimate(X,Y,ymid,F,nsamples,freepar,modelpar,RHS,ode);
else
    % use single error sample per mesh interrval
    nsamples = 1;
    res = residualEstimate(X,Y,ymid,F,nsamples,freepar,modelpar,RHS,ode);
    if max(res) < rtol
        % enter the final verification stage
        OCBVP.solutionVerificationStage = true;
        % compute the remainig two samples
        nsamples = 2;
        res2 = residualEstimate(X,Y,ymid,F,nsamples,freepar,modelpar,RHS,ode);
        res = max(res,res2);
    end
end

function res=residualEstimate(X,Y,ymid,F,nsamples,freepar,modelpar,RHS,ode)
global OCMATCONT OCBVP
sq15=OCBVP.sq15;
nstages=OCBVP.nstages;
threshold=OCBVP.threshold;
switch nsamples
    case 1
        cres = 1/2;
    case 2
        cres = [ 1/2 - sq15/10, 1/2 + sq15/10];
    case 3
        cres = [ 1/2 - sq15/10, 1/2, 1/2 + sq15/10];
end

x = X(1:nstages:end);
y = Y(:,1:nstages:end);
yp = F(:,1:nstages:end);

h = diff(x);

thresh = threshold(:,ones(1,numel(cres)));
res = zeros(size(h));

for i = 1 : numel(h)
    xres = x(i) + cres*h(i);
    [yres,ypres] = ntrp4h(xres,x(i),y(:,i),x(i+1),y(:,i+1),...
        ymid(:,i),yp(:,i),yp(:,i+1));
    residual = ypres - ode(xres,yres,OCMATCONT.HE.arcindex,freepar,modelpar);

    wtdRes = residual ./ max(abs(yres),thresh);
    normRes = max(max(abs(wtdRes)));
    scaledRes = abs(h(i)) * normRes;

    res(i) = scaledRes;
end

