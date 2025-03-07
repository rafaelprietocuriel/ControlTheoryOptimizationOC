function [x,v,Niter] = mynewtcorr2(x0, v0)
%
% [x,v,i] = mynewtcorr2(x0, v0)
% Internal function to perform Newton correction iterations using the
% Newton-Raphson method
% 
% This Newtonalgorithm is based on newtonraphson

global cds
FUN = @(x)funwrapper(x);

x = x0; % initial guess
v = v0;
R = []; 
R(cds.ndim,1) = 1;
%% default values
TYPX = max(abs(x0), 1); % x scaling value, remove zeros
ALPHA = 1e-4; % criteria for decrease
MIN_LAMBDA = 0.1; % min lambda
MAX_LAMBDA = 0.5; % max lambda

%% set scaling values
% TODO: let user set weights
weight = ones(numel(FUN(x0)),1);
J0 = weight*(1./TYPX'); % Jacobian scaling matrix

%% check initial guess
[F,J] = FUN(x); % evaluate initial guess
J=[J;v'];
Jstar = J./J0; % scale Jacobian
if any(isnan(Jstar(:))) || any(isinf(Jstar(:)))
    exitflag = -1; % matrix may be singular
else
    exitflag = 1; % normal exit
end
% if issparse(Jstar)
%     rc = 1/condest(Jstar);
% else
%     if any(isnan(Jstar(:)))
%         rc = NaN;
%     elseif any(isinf(Jstar(:)))
%         rc = Inf;
%     else
%         rc = 1/cond(Jstar); % reciprocal condition
%     end
% end
resnorm = norm(F); % calculate norm of the residuals
dx = zeros(size(x0));
%convergence = Inf; % dummy values
%% solver
Niter = 0; % start counter
lambda = 1; % backtracking
while (resnorm>cds.options.FunTolerance || lambda<1) && exitflag>=0 && Niter<=cds.options.MaxCorrIters
    if lambda==1
        %% Newton-Raphson solver
        Niter = Niter+1; % increment counter
        if cds.options.MoorePenrose % Moore Penrose
            D=Jstar\[F R];
            dx_star=D(:,1);
            v=D(:,2);
            v=v/norm(v);
        else
            dx_star =Jstar\F; % calculate Newton step
        end
        % NOTE: use isnan(f) || isinf(f) instead of STPMAX
        dx = dx_star.*TYPX; % rescale x
        g = F'*Jstar; % gradient of resnorm
        slope = g*dx_star; % slope of gradient
        fold = F'*F; % objective function
        xold = x; % initial value
        lambda_min = cds.options.VarTolerance /max(abs(dx)./max(abs(xold), 1));
    end
    if lambda<lambda_min
        exitflag = 2; % x is too close to XOLD
        break
    elseif any(isnan(dx)) || any(isinf(dx))
        exitflag = -1; % matrix may be singular
        break
    end
    x = xold-dx*lambda; % next guess
    [F, J] = FUN(x); % evaluate next residuals
    J=[J;v'];
    Jstar = J./J0; % scale next Jacobian
    f = F'*F; % new objective function
    %% check for convergence
    lambda1 = lambda; % save previous lambda
    if f>fold+ALPHA*lambda*slope
        if lambda==1
            lambda = -slope/2/(f-fold-slope); % calculate lambda
        else
            A = 1/(lambda1 - lambda2);
            B = [1/lambda1^2,-1/lambda2^2;-lambda2/lambda1^2,lambda1/lambda2^2];
            C = [f-fold-lambda1*slope;f2-fold-lambda2*slope];
            coeff = num2cell(A*B*C);
            [a,b] = coeff{:};
            if a==0
                lambda = -slope/2/b;
            else
                discriminant = b^2 - 3*a*slope;
                if discriminant<0
                    lambda = MAX_LAMBDA*lambda1;
                elseif b<=0
                    lambda = (-b+sqrt(discriminant))/3/a;
                else
                    lambda = -slope/(b+sqrt(discriminant));
                end
            end
            lambda = min(lambda,MAX_LAMBDA*lambda1); % minimum step length
        end
    elseif isnan(f) || isinf(f)
        % limit undefined evaluation or overflow
        lambda = MAX_LAMBDA*lambda1;
    else
        lambda = 1; % fraction of Newton step
    end
    if lambda<1
        lambda2 = lambda1;f2 = f; % save 2nd most previous value
        lambda = max(lambda,MIN_LAMBDA*lambda1); % minimum step length
        continue
    end
    %% display
    %resnorm0 = resnorm; % old resnorm
    %resnorm = norm(F); % calculate new resnorm
    %convergence = log(resnorm0/resnorm); % calculate convergence rate
    %stepnorm = norm(dx); % norm of the step
    if any(isnan(Jstar(:))) || any(isinf(Jstar(:)))
        exitflag = -1; % matrix may be singular
        break
    end
%     if issparse(Jstar)
%         rc = 1/condest(Jstar);
%     else
%         rc = 1/cond(Jstar); % reciprocal condition
%     end
end

if exitflag<0
    x = [];
    v = [];
end
end
function [F, J] = funwrapper(x)
global cds
F=feval(cds.curve_func,x);
F=[F;0];
J=cjac(cds.curve_func,cds.curve_jacobian,x,[]);

end

% BSD 3-Clause License
%
% Copyright (c) 2021, Mark Mikofski
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%
% 1. Redistributions of source code must retain the above copyright notice, this
%    list of conditions and the following disclaimer.
%
% 2. Redistributions in binary form must reproduce the above copyright notice,
%    this list of conditions and the following disclaimer in the documentation
%    and/or other materials provided with the distribution.
%
% 3. Neither the name of the copyright holder nor the names of its
%    contributors may be used to endorse or promote products derived from
%    this software without specific prior written permission.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
% FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
% DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
% SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
% CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
% OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
% OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

