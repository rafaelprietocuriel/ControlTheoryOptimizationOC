function [x,v,iter] = mynewtcorr(x0, v0)
%
% [x,v,i] = mynewtcorr(x0, v0)
% Internal function to perform Newton correction iterations
% 
% This Newtonalgorithm is taken from 'bvp4c'

global cds ics %forpresentation
ics.options.maxnewtiter=4;
ics.options.maxprobes=4;
x = x0;
v = v0;
R = []; R(cds.ndim,1) = 1;
for iter=1:cds.options.MaxCorrIters
    if iter<=cds.options.MaxNewtonIters
        B =cjac(cds.curve_func,cds.curve_jacobian,x,[]);
        wt = max(abs(B),[],2);
        wt(end+1,1)=1;
        B=[B; v'];
        if any(wt == 0) || ~all(isfinite(nonzeros(B)))
            singJac = true;
        else
            scalMatrix = spdiags(1 ./ wt,0,cds.ndim,cds.ndim);
            B = scalMatrix * B;
            [L,U,P] = lu(B);
            singJac=false;
            %singJac = check_singular(dPHIdy,L,U,P,OcSolverVar.warning.warnstat,OcSolverVar.warning.warnoff);
        end
        if singJac
            msg = 'Unable to solve the collocation equations -- a singular Jacobian encountered';
            x=[];
            v=[];
            return
        end
        scalMatrix = P * scalMatrix;
    end
    for t = 1:2
        RHS=[feval(cds.curve_func, x);0];
        if isnan(norm(RHS)) || isinf(norm(RHS))
            x = [];
            v = [];
            return;
        end
        if cds.options.MoorePenrose % Moore Penrose
            lastwarn('');
            % find the Newton direction
            D=U\(L\( scalMatrix * [RHS R] ));
            delxnew=D(:,1);
            vnew=D(:,2);
            v=vnew/norm(vnew);
            if cds.options.CheckSingular &&  ~isempty(lastwarn)
                x = [];
                v = [];
                return;
            end
        else
            delxnew = U\(L\( scalMatrix * RHS ));
        end

        x = x - delxnew;
        if norm(RHS) < cds.options.FunTolerance && norm(delxnew) < cds.options.VarTolerance 
            v = ([cjac(cds.curve_func,cds.curve_jacobian,x,[]);v']\R);
            v = v/norm(v);
            return;
        end
    end
end


x = [];
v = [];

%SD:Newton corrections pal/mp
