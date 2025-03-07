function [xx,yy,tangentnew,NN]=meshadaptation_bvp6c(x,y,tangent,res,canRemovePoints)
%NEW_PROFILE  Redistribute mesh points and approximate the solution.

global OCBVP 

rtol = OCBVP.rtol;
y05=OCBVP.Y05;
Fmid=OCBVP.Fmid;
F=OCBVP.F;
xx = [];
yy = [];
NN = 0;

for arc = 1:OCBVP.numarc

    xidx = OCBVP.Lidx(arc):OCBVP.Ridx(arc);  % mesh point index
    xreg = x(xidx);
    yreg = y(:,xidx);
    Freg = F(:,xidx);
    hreg = diff(xreg);
    Nreg = length(xidx);

    iidx = xidx(1:end-1);    % mesh interval index
    resreg = res(iidx);

    F025reg = Fmid(:,iidx,1);
    F05reg  = Fmid(:,iidx,2);
    F075reg = Fmid(:,iidx,3);

    i1 = find(resreg > rtol);
    i2 = find(resreg(i1) > 100*rtol);
    NNmax = Nreg + length(i1) + length(i2);
    xxreg = zeros(1,NNmax);
    yyreg = zeros(OCBVP.numode,NNmax);
    last_int = Nreg - 1;

    xxreg(1) = xreg(1);
    yyreg(:,1) = yreg(:,1);
    NNreg = 1;
    i = 1;
    while i <= last_int
        if resreg(i) > rtol     % introduce points
            if resreg(i) > 100*rtol
                Ni = OCBVP.numode;
                hi = hreg(i) / (Ni+1);
                j = 1:Ni;
                xxreg(NNreg+j) = xxreg(NNreg) + j*hi;
                for j=1:Ni
                    yyreg(:,NNreg+j) = interp_Hermite_bvp6c(j/(Ni+1),hreg(i),yreg(:,i:i+1),...
                        Freg(:,i:i+1),F025reg(:,i),F05reg(:,i),F075reg(:,i));
                end
            else
                Ni = 1;
                xxreg(NNreg+1) = xxreg(NNreg) + hreg(i)/2;
                yyreg(:,NNreg+1) = y05(:,i);
            end
            NNreg = NNreg + Ni;
        else                 % try removing points
            if canRemovePoints && (i <= last_int-2) && (max(resreg(i+1:i+2)) < rtol)
                hnew = (hreg(i)+hreg(i+1)+hreg(i+2))/2;
                C1 = resreg(i)/(hreg(i)/hnew)^(11/2);
                C2 = resreg(i+1)/(hreg(i+1)/hnew)^(11/2);
                C3 = resreg(i+2)/(hreg(i+2)/hnew)^(11/2);
                pred_res = max([C1,C2,C3]);

                if pred_res < 0.05 * rtol   % replace 3 intervals with 2
                    xxreg(NNreg+1) = xxreg(NNreg) + hnew;
                    yyreg(:,NNreg+1) = interp_Hermite_bvp6c(0.5,hreg(i),yreg(:,(i+1):(i+2)),...
                        Freg(:,(i+1):(i+2)),F025reg(:,i+1),F05reg(:,i+1),F075reg(:,i+1));
                    NNreg = NNreg + 1;
                    i = i + 2;
                end
            end
        end
        NNreg = NNreg + 1;
        xxreg(NNreg) = xreg(i+1);   % preserve the next mesh point
        yyreg(:,NNreg) = yreg(:,i+1);
        i = i + 1;
    end

    NN = NN + NNreg;
    if (NN > OCBVP.Nmax)
        % return the previous solution
        xx = x;
        yy = y;
        %mbcidxnew = mbcidx;
        break
    else
        xx = [xx, xxreg(1:NNreg)];
        yy = [yy, yyreg(:,1:NNreg)];
        %if region < nregions    % possible only for multipoint BVPs
        %    mbcidxnew = [mbcidxnew, NN];
        %end
    end
end

Nold=length(x);
nNold=OCBVP.numode*Nold;
tangentnew=devaltangent(xx,x,tangent,OCBVP.numode,Nold,nNold);
