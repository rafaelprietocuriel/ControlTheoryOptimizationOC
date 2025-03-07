function [coeff,tangent,iter]=semismoothnewtcorr4stat(coeff0,tangent0)
% Newton solver and mesh adaptation

global OCSTATCONT

needGlobJac=true;

coeff=coeff0;
tangent=tangent0(:).';
EN(OCSTATCONT.HE.numdvariables,1)=1;

% THE MAIN LOOP:
iter=0;
PHI=OCSTATCONT.generalizedoperatoreq(coeff,tangent);
for iter=1:OCSTATCONT.OPTIONS.maxnewtiter
    if needGlobJac
        DPHI=OCSTATCONT.generalizedfrechetder(coeff,tangent);
        DPHI(OCSTATCONT.HE.numdvariables,:)=tangent(:).';
    end
    if OCSTATCONT.OPTIONS.continuationmethod==1 % Moore Penrose
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
    lambda = 1;
    for probe = 1:OCSTATCONT.OPTIONS.maxprobes
        coeff_new = coeff-lambda*coeff_del;
        PHI=OCSTATCONT.generalizedoperatoreq(coeff_new,tangent);
        coeff_new_dist=norm(DPHI\PHI);

        if (coeff_new_dist < 0.9*coeff_dist)
            break
        else
            lambda = 0.5*lambda;
        end
    end
    needGlobJac = (coeff_new_dist > 0.1*coeff_dist);
    coeff = coeff_new;
    if coeff_new_dist < 0.1*OCSTATCONT.OPTIONS.newtonabstol
        DPHI=OCSTATCONT.generalizedfrechetder(coeff,tangent);
        DPHI(OCSTATCONT.HE.numdvariables,:)=tangent(:).';
        tangent=DPHI\EN;
        tangent=tangent/norm(tangent);
        OCSTATCONT.RHS=PHI;
        return
    end
end
coeff=[];
tangent=[];
OCSTATCONT.RHS=[];
OCSTATCONT.J=[];