function maxInterpResidual=maxresidual_bvp5c(X,Y)
global OCBVP
%

% This function bases on bvp4c by
%   Jacek Kierzenka and Lawrence F. Shampine
%   Copyright 1984-2007 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $  $Date: 2007/05/23 18:54:08 $

F=OCBVP.F;
B=OCBVP.B;
sq5=OCBVP.sq5;
nstages=OCBVP.nstages;
threshold=OCBVP.threshold;
maxInterpResidual = 0;

h = diff(X(1:nstages:end));
thresh = threshold(:,ones(1,nstages));

idx = 0;
for i = 1 : numel(h)
    ync = Y(:,  idx+2 : idx+nstages+1);
    K   = F(:,  idx+1 : idx+nstages+1);

    % Form the residual in the collocation equations in terms of the
    % intermediate derivative values and compute a weighted norm scaled
    % by the mesh spacing.
    K1 = K(:,1);
    Boh = B/h(i);
    K2_4 = [ -sq5/5*K1, sq5/5*K1, -K1] + Y(:,idx+1:idx+4)*Boh';

    residual = K2_4 - K(:,2:4);

    wtdRes = residual ./ max(abs(ync),thresh);
    normRes = max(max(abs(wtdRes)));
    scaledRes = abs(h(i)) * normRes;

    maxInterpResidual = max(maxInterpResidual,scaledRes);

    idx = idx + nstages;
end
