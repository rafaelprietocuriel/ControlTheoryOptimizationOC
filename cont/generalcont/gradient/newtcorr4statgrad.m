function [coeff,extremal,tangent,iter,maxgraditer]=newtcorr4statgrad(coeff0,extremal0,tangent0)
% Newton solver and mesh adaptation

global OCGRADCONT
coeff=coeff0;
tangent=tangent0(:).';
extremal=extremal0;
maxgraditer=-inf;
EN(OCGRADCONT.HE.numdvariables,1)=1;

% THE MAIN LOOP:
for iter=1:OCGRADCONT.OPTIONS.maxnewtiter
    if iter<=OCGRADCONT.OPTIONS.maxnewtiter
        DPHI=OCGRADCONT.frechetder(coeff,extremal,tangent);
        if isempty(DPHI)
            [extremal,maxgraditer]=OCGRADCONT.gradientsolution(coeff,extremal,tangent);
            iter=0;
            return
        end
        DPHI(OCGRADCONT.HE.numdvariables,:)=tangent(:).';
    end
    for probe = 1:OCGRADCONT.OPTIONS.maxprobes
        [PHI,extremal,graditer]=OCGRADCONT.operatoreq(coeff,extremal,tangent);
        maxgraditer=max(maxgraditer,graditer);
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
        coeff = coeff - coeff_del;
        if norm(coeff_del) < OCGRADCONT.OPTIONS.newtonreltol && norm(PHI) < OCGRADCONT.OPTIONS.newtonabstol
            DPHI=OCGRADCONT.frechetder(coeff,extremal,tangent);
            DPHI(OCGRADCONT.HE.numdvariables,:)=tangent(:).';
            tangent=DPHI\EN;
            tangent=tangent/norm(tangent);
            OCGRADCONT.RHS=PHI;
            return
        end
    end
end
coeff=[];
tangent=[];
OCGRADCONT.RHS=[];
OCGRADCONT.J=[];