function [xx,yy,tangentnew,NN]=meshadaptation_gbvp4c(x,y,tangent,res,canRemovePoints)
%NEW_PROFILE  Redistribute mesh points and approximate the solution.

% This function bases on bvp4c by
%   Jacek Kierzenka and Lawrence F. Shampine
%   Copyright 1984-2007 The MathWorks, Inc.
%   $Revision: 1.21.4.15 $  $Date: 2007/05/23 18:54:07 $

global OCBVP 

yp=OCBVP.F;
rtol = OCBVP.rtol;

xx = [];
yy = [];
NN = 0;

for arc = 1:OCBVP.numarc

    xidx = OCBVP.Lidx(arc):OCBVP.Ridx(arc);  % mesh point index
    xreg = x(xidx);
    yreg = y(1:OCBVP.numode(arc),xidx);
    ypreg = yp(1:OCBVP.numode(arc),xidx);
    hreg = diff(xreg);
    Nreg = length(xidx);

    iidx = xidx(1:end-1);    % mesh interval index
    resreg = res(iidx);
    i1 = find(resreg > rtol);
    i2 = find(resreg(i1) > 100*rtol);
    NNmax = Nreg + length(i1) + length(i2);
    xxreg = zeros(1,NNmax);
    yyreg = zeros(OCBVP.numode(arc),NNmax);
    last_int = Nreg - 1;

    xxreg(1) = xreg(1);
    yyreg(:,1) = yreg(:,1);
    NNreg = 1;
    i = 1;
    while i <= last_int
        if resreg(i) > rtol     % introduce points
            if resreg(i) > 100*rtol
                Ni = 2;
            else
                Ni = 1;
            end
            hi = hreg(i) / (Ni+1);
            j = 1:Ni;
            xxreg(NNreg+j) = xxreg(NNreg) + j*hi;
            yyreg(:,NNreg+j) = ntrp3h(xxreg(NNreg+j),xreg(i),yreg(:,i),xreg(i+1),...
                yreg(:,i+1),ypreg(:,i),ypreg(:,i+1));
            NNreg = NNreg + Ni;
        else
            if canRemovePoints && (i <= last_int-2) && all(resreg(i+1:i+2) < rtol)
                % try removing points
                hnew = (hreg(i)+hreg(i+1)+hreg(i+2))/2;
                C1 = resreg(i)/(hreg(i)/hnew)^(7/2);
                C2 = resreg(i+1)/(hreg(i+1)/hnew)^(7/2);
                C3 = resreg(i+2)/(hreg(i+2)/hnew)^(7/2);
                pred_res = max([C1,C2,C3]);

                if pred_res < 0.5 * rtol   % replace 3 intervals with 2
                    xxreg(NNreg+1) = xxreg(NNreg) + hnew;
                    yyreg(:,NNreg+1) = ntrp3h(xxreg(NNreg+1),xreg(i+1),yreg(:,i+1),xreg(i+2),...
                        yreg(:,i+2),ypreg(:,i+1),ypreg(:,i+2));
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
        yyregnew=zeros(OCBVP.maxnumode,NNreg);
        yyregnew(1:OCBVP.numode(arc),:)=yyreg(:,1:NNreg);
        yy = [yy, yyregnew];
        %if region < nregions    % possible only for multipoint BVPs
        %    mbcidxnew = [mbcidxnew, NN];
        %end
    end
end
Nold=length(x);
nNold=OCBVP.maxnumode*Nold;
tangentnew=devaltangent(xx,x,tangent,OCBVP.numode,Nold,nNold);

function [yint,ypint] = ntrp3h(tint,t,y,tnew,ynew,yp,ypnew)
%NTRP3H  Interpolation helper function for BVP4C, DDE23, and DDESD.
%   YINT = NTRP3H(TINT,T,Y,TNEW,YNEW,YP,YPNEW) evaluates the Hermite cubic
%   interpolant at time TINT. TINT may be a scalar or a row vector.
%   [YINT,YPINT] = NTRP3H(TINT,T,Y,TNEW,YNEW,YP,YPNEW) returns also the
%   derivative of the interpolating polynomial.
%
%   See also BVP4C, DDE23, DDESD, DEVAL.

%   Jacek Kierzenka and Lawrence F. Shampine
%   Copyright 1984-2005 The MathWorks, Inc.
%   $Revision: 1.1.6.5 $  $Date: 2005/11/18 14:15:50 $

h = tnew - t;
s = (tint - t)/h;
s2 = s .* s;
s3 = s .* s2;
slope = (ynew - y)/h;
c = 3*slope - 2*yp - ypnew;
d = yp + ypnew - 2*slope;
yint = y(:,ones(size(tint))) + (h*d*s3 + h*c*s2 + h*yp*s);
if nargout > 1
    ypint = yp(:,ones(size(tint))) + (3*d*s2 + 2*c*s);
end

