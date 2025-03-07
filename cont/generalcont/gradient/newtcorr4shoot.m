function [coeff,extremal,tangent,iter]=newtcorr4shoot(coeff0,extremal0,tangent0)
% Newton solver and mesh adaptation

global OCSHOOTCONT
needGlobJac=true;
coeff=coeff0;
tangent=tangent0(:).';
extremal=extremal0;
EN=OCSHOOTCONT.HE.unitvector;
% THE MAIN LOOP:
[PHI,extremal]=OCSHOOTCONT.operatoreq(coeff,extremal,tangent);
for iter=1:OCSHOOTCONT.OPTIONS.maxnewtiter
    if needGlobJac
        % setup and factor the global Jacobian
        DPHI=OCSHOOTCONT.frechetder(coeff,extremal,tangent);
        DPHI(OCSHOOTCONT.HE.numdvariables,:)=tangent(:).';
        %needGlobJac = false;
    end
    if OCSHOOTCONT.OPTIONS.continuationmethod==1 % Moore Penrose
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
    for probe = 1:OCSHOOTCONT.OPTIONS.maxprobes
        coeff_new = coeff - lambda*coeff_del;
        extremal_new=extremal;%OCSHOOTCONT.predictextremal(coeff_new,extremal,-lambda*coeff);
        [PHI,extremal_new]=OCSHOOTCONT.operatoreq(coeff_new,extremal_new,tangent);

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
    if coeff_new_dist < 0.1*OCSHOOTCONT.OPTIONS.newtonabstol
        DPHI=OCSHOOTCONT.frechetder(coeff,extremal,tangent);
        DPHI(OCSHOOTCONT.HE.numdvariables,:)=tangent(:).';
        tangent=DPHI\EN;
        tangent=tangent/norm(tangent);
        OCSHOOTCONT.RHS=PHI;
        return
    end
end
coeff=[];
tangent=[];
OCSHOOTCONT.RHS=[];
OCSHOOTCONT.J=[];