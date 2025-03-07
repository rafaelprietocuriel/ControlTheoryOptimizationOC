function [tmesh,coeff,iter]=newtcorr4bvpsolve(tmesh0,coeff0)
% Newton solver and mesh adaptation

global OCMATCONT OCBVP
OCBVP.Joptions.fac=[];
OCBVP.dPoptions.fac=[];
needGlobJac=true;
tangent=[];
coeff=coeff0;
tmesh=tmesh0;


% THE MAIN LOOP:
PHI=OCMATCONT.operatoreq(tmesh,coeff,tangent,OCMATCONT.ode,OCMATCONT.bc,OCMATCONT.icfun);
for iter=1:OCMATCONT.OPTIONS.maxnewtiter
    if needGlobJac
        % setup and factor the global Jacobian
        DPHI=OCMATCONT.frechetder(tmesh,coeff,tangent,OCMATCONT.ode,OCMATCONT.bc,OCMATCONT.icfun,OCMATCONT.odejac,OCMATCONT.bcjac,OCMATCONT.icfunjac);
        wt = max(abs(DPHI),[],2);
        needGlobJac = true;
        if any(wt == 0) || ~all(isfinite(nonzeros(DPHI)))
            singJac = true;
        else
            %scalMatrix =speye(OCMATCONT.HE.numdvariables);
            scalMatrix = spdiags(1 ./ wt,0,OCMATCONT.HE.numdvariables,OCMATCONT.HE.numdvariables);
            DPHI = scalMatrix * DPHI;
            [L,U,P] = lu(DPHI);
            if OCMATCONT.OPTIONS.checksingular
                singJac = check_singular(DPHI,L,U,P);
            else
                singJac=false;
            end
        end
        if singJac
            ocmatmsg('Unable to solve the collocation equations -- a singular Jacobian encountered\n');
            coeff=[];
            return
        end
        scalMatrix = P * scalMatrix;
    end
    coeff_del = U\(L\(scalMatrix*PHI));
    coeff_dist = norm(coeff_del);
    % weak line search with an affine stopping criterion
    lambda = 1;
    for probe = 1:OCMATCONT.OPTIONS.maxprobes
        coeff_new = coeff - lambda*coeff_del;
        PHI=OCMATCONT.operatoreq(tmesh,coeff_new,tangent,OCMATCONT.ode,OCMATCONT.bc,OCMATCONT.icfun);
        OCMATCONT.RHS=PHI;
        coeff_new_dist = norm(U\(L\( scalMatrix * PHI )));
        %coeff_new_dist = norm(U\(L\( PHI )));

        if (coeff_new_dist < 0.9*coeff_dist)
            break
        else
            lambda = 0.5*lambda;
        end
    end

    needGlobJac = (coeff_new_dist > 0.1*coeff_dist);
    coeff = coeff_new;
    if coeff_new_dist < OCMATCONT.OPTIONS.newtonabstol
        return
    end
end
tmesh=[];
coeff=[];
