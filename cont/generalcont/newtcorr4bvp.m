function [tmesh,coeff,tangent,iter]=newtcorr4bvp(tmesh0,coeff0,tangent0)
% Newton solver and mesh adaptation

global OCMATCONT OCBVP
OCBVP.Joptions.fac=[];
OCBVP.dPoptions.fac=[];
needGlobJac=true;

coeff=coeff0;
tangent=tangent0(:).';
tmesh=tmesh0;

EN=[];
EN(OCMATCONT.HE.numdvariables,1)=1;

wt=zeros(OCMATCONT.HE.numdvariables,1);
% THE MAIN LOOP:
PHI=OCMATCONT.operatoreq(tmesh,coeff,tangent,OCMATCONT.ode,OCMATCONT.bc,OCMATCONT.icfun);
for iter=1:OCMATCONT.OPTIONS.maxnewtiter
    if needGlobJac
        % setup and factor the global Jacobian
        DPHI=OCMATCONT.frechetder(tmesh,coeff,tangent,OCMATCONT.ode,OCMATCONT.bc,OCMATCONT.icfun,OCMATCONT.odejac,OCMATCONT.bcjac,OCMATCONT.icfunjac);
        if isempty(DPHI)
            break
        end
        wt(1:OCMATCONT.HE.numdvariables-1) = max(abs(DPHI),[],2);
        % last row is not rescaled since it consists of the orthogonal
        % search direction
        wt(OCMATCONT.HE.numdvariables,1)=1;
        DPHI=[DPHI;tangent(:).'];
        %needGlobJac = false;
        if any(wt == 0) || ~all(isfinite(nonzeros(DPHI)))
            singJac = true;
        else
            scalMatrix = spdiags(1 ./ wt,0,OCMATCONT.HE.numdvariables,OCMATCONT.HE.numdvariables);
            DPHI = scalMatrix * DPHI;
            %clear L U P
            [L,U,P] = lu(DPHI);
            if  OCMATCONT.OPTIONS.checksingular
                singJac = check_singular(DPHI,L,U,P);
            else
                singJac=false;
            end
        end
        if singJac
            ocmatmsg('Unable to solve the collocation equations -- a singular Jacobian encountered\n');
            coeff=[];
            tangent=[];
            return
        end
        scalMatrix = P * scalMatrix;
    end
    if OCMATCONT.OPTIONS.continuationmethod==1 % Moore Penrose
        % find the Newton direction
        D=U\(L\(scalMatrix*[PHI EN] ));
        coeff_del=D(:,1);
        tangent=D(:,2);
        tangent=tangent/norm(tangent);
    else % Arclength
        coeff_del = U\(L\(scalMatrix*PHI));
    end
    coeff_dist = norm(coeff_del);
    % weak line search with an affine stopping criterion
    lambda = 1;
    for probe = 1:OCMATCONT.OPTIONS.maxprobes
        coeff_new = coeff - lambda*coeff_del;
        PHI=OCMATCONT.operatoreq(tmesh,coeff_new,tangent,OCMATCONT.ode,OCMATCONT.bc,OCMATCONT.icfun);
        coeff_new_dist = norm(U\(L\( scalMatrix * PHI )));

        if (coeff_new_dist < 0.9*coeff_dist)
            break
        else
            lambda = 0.5*lambda;
        end
    end

    needGlobJac = (coeff_new_dist > 0.1*coeff_dist);

    coeff = coeff_new;
    if coeff_new_dist<OCMATCONT.OPTIONS.newtonreltol
        DPHI=OCMATCONT.frechetder(tmesh,coeff,tangent,OCMATCONT.ode,OCMATCONT.bc,OCMATCONT.icfun,OCMATCONT.odejac,OCMATCONT.bcjac,OCMATCONT.icfunjac);
        wt = max(abs(DPHI),[],2);
        % last row is not rescaled since it consists of the orthogonal
        % search direction
        wt(end+1,1)=1;
        scalMatrix = spdiags(1 ./ wt,0,OCMATCONT.HE.numdvariables,OCMATCONT.HE.numdvariables);
        DPHI=[DPHI;tangent(:).'];
        DPHI = scalMatrix * DPHI;
        [L,U,P] = lu(DPHI);
        scalMatrix = P * scalMatrix;
        tangent=U\(L\(scalMatrix*EN));
        %         scalMatrix = spdiags(1 ./ wt,0,OCMATCONT.HE.numdvariables,OCMATCONT.HE.numdvariables);
        %         DPHI = scalMatrix * DPHI;
        %         clear L U P
        %         [L,U,P] = lu(DPHI);
        %         scalMatrix = P * scalMatrix;
        %         tangent=U\(L\(scalMatrix*EN));
        tangent=tangent/norm(tangent);
        OCMATCONT.RHS=PHI;
        %OCMATCONT.J=DPHI;
        return
    end
end
tmesh=[];
coeff=[];
tangent=[];
OCMATCONT.RHS=[];
OCMATCONT.J=[];
