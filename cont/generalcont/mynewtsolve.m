function [x,fval,exitflag,info] = mynewtsolve(x0,zerofunc,jacfunc,opt,varargin)
%
% [x,v,i] = mynewtcorr(x0, v0)
% Internal function to perform Newton correction iterations
%
% This Newtonalgorithm is taken from 'bvp4c'

if isempty(opt)
    opt=defaultocoptions;
end
x = x0;
dim=length(x);
for iter=1:opt.NEWTON.MaxProbes
    if iter<=opt.NEWTON.MaxNewtonIters
        B=calcjacobian(x,varargin{:});
        wt = max(abs(B),[],2);
        if any(wt == 0) || ~all(isfinite(nonzeros(B)))
            singJac = true;
        else
            scalMatrix = spdiags(1 ./ wt,0,dim,dim);
            B = scalMatrix * B;
            [L,U,P] = lu(B);
            singJac=false;
            %singJac = check_singular(dPHIdy,L,U,P,OcSolverVar.warning.warnstat,OcSolverVar.warning.warnoff);
        end
        if singJac
            ocmatmsg('Unable to solve the collocation equations -- a singular Jacobian encountered\n');
            exitflag=0;
            fval=zerofunc(x,varargin{:});
            info.Iterations=iter;
            return
        end
        scalMatrix = P * scalMatrix;
    end
    for t = 1:2
        RHS=zerofunc(x,varargin{:});
        if isnan(norm(RHS)) || isinf(norm(RHS))
            x = [];
            exitflag=0;
            fval=[];
            info.Iterations=iter;
            return;
        end
        delxnew = U\(L\( scalMatrix * RHS ));

        x = x - delxnew;
        if norm(RHS) < opt.NEWTON.AbsTol && norm(delxnew) < opt.NEWTON.RelTol
            fval=zerofunc(x,varargin{:});
            exitflag=1;
            info.Iterations=iter;
            return;
        end
    end
end
exitflag=0;
fval=zerofunc(x,varargin{:});
info.Iterations=iter;

    function J=calcjacobian(x,varargin)
        if isempty(jacfunc)
            numJacOpt.diffvar=1;
            numJacOpt.vectvars=[];
            J=numjaccsd(zerofunc,{x,varargin{:}},length(x),numJacOpt);
        else
            J=jacfunc(x,varargin{:});
        end
    end
end