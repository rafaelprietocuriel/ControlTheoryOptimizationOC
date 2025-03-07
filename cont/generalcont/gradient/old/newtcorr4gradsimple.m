function [coeff,extremal,tangent,iter,maxgraditer,maxdgraditer,maxlineiter,maxdlineiter,numrepdd]=newtcorr4gradsimple(coeff0,extremal0,tangent0)
% Newton solver and mesh adaptation

global OCGRADCONT
needGlobJac=true;
coeff=coeff0;
tangent=tangent0(:).';
extremal=extremal0;
maxgraditer=-inf;
maxdgraditer=-inf;
maxlineiter=-inf;
maxdlineiter=-inf;
EN=OCGRADCONT.HE.unitvector;
% THE MAIN LOOP:
[PHI,extremal,graditer,dgraditer,lineiter,dlineiter,numrepdd]=OCGRADCONT.operatoreq(coeff,extremal,tangent);
if isempty(PHI)
    [extremal,maxgraditer,maxdgraditer,maxlineiter,maxdlineiter,numrepdd]=OCGRADCONT.gradientsolution(coeff,extremal,tangent);
    iter=0;
    return
end
for iter=1:OCGRADCONT.OPTIONS.maxnewtiter
    if needGlobJac
        % setup and factor the global Jacobian
        DPHI=OCGRADCONT.frechetder(coeff,extremal,tangent);
        DPHI(OCGRADCONT.HE.numdvariables,:)=tangent(:).';
        %needGlobJac = false;
    end
    if OCGRADCONT.OPTIONS.continuationmethod==1 % Moore Penrose
        % find the Newton direction
        D=DPHI\[PHI EN];
        lastwarn('');
        if ~isempty(lastwarn)
            coeff=[];
            tangent=[];
            return;
        end
        coeff_del=D(:,1);
        tangent=D(:,2);
        tangent=tangent/norm(tangent);
    else % Arclength
        coeff_del = DPHI\PHI;
    end
    coeff_dist = norm(coeff_del);
    % weak line search with an affine stopping criterion
    lambda = 1;
    for probe = 1:OCGRADCONT.OPTIONS.maxprobes
        coeff_new = coeff - lambda*coeff_del;
        extremal_new=extremal;
        [PHI,extremal_new,graditer,dgraditer,lineiter,dlineiter,numrepdd]=OCGRADCONT.operatoreq(coeff_new,extremal_new,tangent);
        maxgraditer=max(maxgraditer,graditer);
        maxdgraditer=max(maxdgraditer,dgraditer);
        maxlineiter=max(maxlineiter,lineiter);
        maxdlineiter=max(maxdlineiter,dlineiter);

        coeff_new_dist = norm(PHI);

        if (coeff_new_dist < 0.9*coeff_dist)
            break
        else
            lambda = 0.5*lambda;
        end
    end

    needGlobJac = (coeff_new_dist > 0.1*coeff_dist);

    coeff = coeff_new;
    extremal=extremal_new;
    if coeff_new_dist < 0.1*OCGRADCONT.OPTIONS.newtonabstol
        DPHI=OCGRADCONT.frechetder(coeff,extremal,tangent);
        DPHI(OCGRADCONT.HE.numdvariables,:)=tangent(:).';
        tangent=DPHI\EN;
        tangent=tangent/norm(tangent);
        OCGRADCONT.RHS=PHI;
        return
    end
end
coeff=[];
tangent=[];
OCGRADCONT.RHS=[];
OCGRADCONT.J=[];