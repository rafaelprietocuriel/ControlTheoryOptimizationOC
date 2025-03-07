function [XX,YY,tangentnew,NN]=meshadaptation_bvp5c(X,Y,tangent,res,canRemovePoints,freepar,modelpar,ode)
%NEW_PROFILE  Redistribute mesh points and approximate the solution.

global OCMATCONT OCBVP

F=OCBVP.F;
nstages=OCBVP.nstages;
rtol=OCBVP.rtol;
ymid=OCBVP.ymid;
c=OCBVP.c;
neqn=OCBVP.neqn;

pow = 5;

% extract the derivative at mesh points
yp = F(:,1:nstages:end);

XX(1) = X(1);
YY(:,1) = Y(:,1);
FF(:,1) = F(:,1);

dirx = sign(X(end)-X(1));
xidx = 1;
i = 1;

while i <= numel(res)

    if res(i) <= rtol

        if canRemovePoints && ...
                ((i+2) <= numel(res)) && all(res(i:i+2) < 0.1*rtol)
            % consider consolidating mesh intervals

            xi = X(xidx);
            yi = Y(:,xidx);
            xip1 = X(xidx+nstages);
            yip1 = Y(:,xidx+nstages);
            xip2 = X(xidx+2*nstages);
            yip2 = Y(:,xidx+2*nstages);
            xip3 = X(xidx+3*nstages);
            yip3 = Y(:,xidx+3*nstages);
            fip3 = F(:,xidx+3*nstages);
            hnew = (xip3 - xi)/2;

            % predict new residuals
            C1 = res(i)   / abs(xip1-xi)^pow;
            C2 = res(i+1) / abs(xip2-xip1)^pow;
            C3 = res(i+2) / abs(xip3-xip2)^pow;
            predEst = max([C1,C2,C3])*abs(hnew)^pow;

            if predEst < 0.5 * rtol  % replace 3 intervals with 2
                xx = xi + hnew*[c,1,(1+c)];
                yy = zeros(neqn,numel(xx));
                % interpolate xx from xi,xip1; then from xip1,xip2; then from xip2,xip3
                idx = find(dirx*(xx - xip1) <= 0);
                if ~isempty(idx)
                    yy(:,idx) = ntrp4h(xx(idx),xi,yi,xip1,yip1,ymid(:,i),yp(:,i),yp(:,i+1));
                end
                idx = find( (dirx*(xx -xip1) > 0) & (dirx*(xx - xip2) <= 0));
                if ~isempty(idx)
                    yy(:,idx) = ntrp4h(xx(idx),xip1,yip1,xip2,yip2,ymid(:,i+1),yp(:,i+1),yp(:,i+2));
                end
                idx = find( dirx*(xx - xip2) > 0);
                if ~isempty(idx)
                    yy(:,idx) = ntrp4h(xx(idx),xip2,yip2,xip3,yip3,ymid(:,i+2),yp(:,i+2),yp(:,i+3));
                end
                ff = ode(xx,yy,OCMATCONT.HE.arcindex,freepar,modelpar);

                xx(  end+1) = xip3;
                yy(:,end+1) = yip3;
                ff(:,end+1) = fip3;

                i = i + 3;
                xidx = xidx + 3*nstages;

            else
                xx = X(   xidx+1 : xidx+nstages);
                yy = Y(:, xidx+1 : xidx+nstages);
                ff = F(:, xidx+1 : xidx+nstages);
                i = i + 1;
                xidx = xidx + nstages;

            end

        else
            xx = X(   xidx+1 : xidx+nstages);
            yy = Y(:, xidx+1 : xidx+nstages);
            ff = F(:, xidx+1 : xidx+nstages);
            i = i + 1;
            xidx = xidx + nstages;

        end

    else  % res(i) > rtol
        % split mesh interval
        xi = X(xidx);
        yi = Y(:,xidx);
        xip1 = X(  xidx+nstages);
        yip1 = Y(:,xidx+nstages);
        fip1 = F(:,xidx+nstages);

        if res(i) > 250 * rtol
            % split into three -- introduce two points
            hnew = (xip1 - xi)/3;
            xx = xi + hnew*[c,1,(1+c),2,(2+c)];
        else
            % split into two -- introduce one point
            hnew = (xip1 - xi)/2;
            xx = xi + hnew*[c,1,(1+c)];
        end
        yy = ntrp4h(xx,xi,yi,xip1,yip1,ymid(:,i),yp(:,i),yp(:,i+1));
        ff =  ode(xx,yy,OCMATCONT.HE.arcindex,freepar,modelpar);
        xx(  end+1) = xip1;
        yy(:,end+1) = yip1;
        ff(:,end+1) = fip1;
        i = i + 1;
        xidx = xidx + nstages;
    end
    XX(   end+1 : end+numel(xx)) = xx;
    YY(:, end+1 : end+numel(xx)) = yy;
    FF(:, end+1 : end+numel(xx)) = ff;
end

OCBVP.F=FF;
NN=numel(XX);
Nold=length(x);
nNold=OCBVP.numode*Nold;
tangentnew=devaltangent(xx,x,tangent,OCBVP.numode,Nold,nNold);

