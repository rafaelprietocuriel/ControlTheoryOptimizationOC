function [tmesh,coeff,tangent,iter]=newtcorr4bvpI(tmesh0,coeff0,tangent0)
% Newton solver and mesh adaptation

global OCMATCONT 
needGlobJac=true;

coeff=coeff0;
tangent=tangent0(:).';
tmesh=tmesh0;

meshHistory = [0,0];            % Keep track of [N, maxres],

EN=[];
EN(OCMATCONT.HE.numdvariables,1)=1;
done=false;
% THE MAIN LOOP:
while ~done
    PHI=OCMATCONT.operatoreq(tmesh,coeff,tangent,OCMATCONT.ode,OCMATCONT.bc,OCMATCONT.icfun);
    for iter=1:OCMATCONT.OPTIONS.maxnewtiter
        if needGlobJac
            % setup and factor the global Jacobian
            DPHI=OCMATCONT.frechetder(tmesh,coeff,tangent,OCMATCONT.ode,OCMATCONT.bc,OCMATCONT.icfun,OCMATCONT.odejac,OCMATCONT.bcjac,OCMATCONT.icfunjac);
            wt = max(abs(DPHI),[],2);
            % last row is not rescaled it consists of the orthogonal search
            % direction
            wt(end+1,1)=1;
            DPHI=[DPHI;tangent(:).'];
            needGlobJac = false;
            if any(wt == 0) || ~all(isfinite(nonzeros(DPHI)))
                singJac = true;
            else
                scalMatrix = spdiags(1 ./ wt,0,OCMATCONT.HE.numdvariables,OCMATCONT.HE.numdvariables);
                %scalMatrix =speye(OCMATCONT.HE.numdvariables);
                DPHI = scalMatrix * DPHI;
                [L,U,P,Q,R] = lu(DPHI);
                singJac = bvpcheckmatrix(DPHI,L,U,P,Q,R);
            end
            if singJac
                ocmatmsg('Unable to solve the collocation equations -- a singular Jacobian encountered\n');
                coeff=[];
                tangent=[];
                return
            end
        end
        if OCMATCONT.OPTIONS.continuationmethod==1 % Moore Penrose
            % find the Newton direction
            D=Q * (U \ (L \ (P * (R \ ( scalMatrix*[PHI EN])))));
            coeff_del=D(:,1);
            tangent=D(:,2);
            tangent=tangent/norm(tangent);
        else % Arclength
            coeff_del = Q * (U \ (L \ (P * (R \ (scalMatrix*PHI)))));
        end
        coeff_dist = norm(coeff_del);
        % weak line search with an affine stopping criterion
        lambda = 1;
        for probe = 1:OCMATCONT.OPTIONS.maxprobes
            coeff_new = coeff - lambda*coeff_del;
            PHI=OCMATCONT.operatoreq(tmesh,coeff_new,tangent,OCMATCONT.ode,OCMATCONT.bc);
            coeff_new_dist = norm(U \ (L \ (P * (R \ (scalMatrix * PHI)))));

            if (coeff_new_dist < 0.9*coeff_dist)
                break
            else
                lambda = 0.5*lambda;
            end
        end

        needGlobJac = (coeff_new_dist > 0.1*coeff_dist);

        coeff = coeff_new;
        if coeff_new_dist < 0.1*OCMATCONT.OPTIONS.newtonabstol
            if OCMATCONT.maxresidual(tmesh,coeff)< 0.1*OCMATCONT.OPTIONS.newtonabstol
                DPHI=OCMATCONT.frechetder(tmesh,coeff,tangent,OCMATCONT.ode,OCMATCONT.bc,OCMATCONT.odejac,OCMATCONT.bcjac);
                DPHI=[DPHI;tangent(:).'];
                [L,U] = lu(DPHI);
                tangent=U\(L\EN);
                tangent=tangent/norm(tangent);
                break
            end
        end
    end
    [tmesh,coeff,tangent,done,meshHistory,needGlobJac,EN]= ...
        meshadaptation(tmesh,coeff,tangent,PHI,done,meshHistory,needGlobJac,EN);
end
