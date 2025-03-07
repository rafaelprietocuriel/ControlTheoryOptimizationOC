function [F Fmid]=calcdata_bvp4c(x,y,freepar,modelpar,ode)

global OCMATCONT OCBVP

% multi-point BVP support
arcposition = find(diff(x) == 0);
Lidx = [1, arcposition+1];
Ridx = [arcposition, length(x)];
FcnArgs = {0,freepar,modelpar};   
n=OCBVP.numode;
N=numel(x);

F = zeros(n,N);
Fmid = zeros(n,N-1);    % include interface intervals

for arc = 1:OCBVP.numarc
    FcnArgs{1} = arc;
    xidx = Lidx(arc):Ridx(arc);   % mesh point index
    Nreg = length(xidx);
    xreg = x(xidx);
    yreg = y(:,xidx);

    iidx = xidx(1:end-1);   % mesh interval index
    Nint = length(iidx);
    % derivative at the mesh points
    if OCMATCONT.OPTIONS.xyvectorized
        Freg = ode(xreg,yreg,FcnArgs{:});
    else
        Freg = zeros(n,Nreg);
        for i = 1:Nreg
            Freg(:,i) = ode(xreg(i),yreg(:,i),FcnArgs{:});
        end
    end
    % derivative at the midpoints
    h = diff(xreg);
    H = spdiags(h(:),0,Nint,Nint);
    xi = xreg(1:end-1);
    yi = yreg(:,1:end-1);
    xip1 = xreg(2:end);
    yip1 = yreg(:,2:end);
    Fi = Freg(:,1:end-1);
    Fip1 = Freg(:,2:end);

    xip05 = (xi + xip1)/2;
    yip05 = (yi + yip1)/2 - (Fip1 - Fi)*(H/8);
    if OCMATCONT.OPTIONS.xyvectorized
        Fip05 = ode(xip05,yip05,FcnArgs{:});
    else % not vectorized
        Fip05 = zeros(n,Nint);
        for i = 1:Nint
            Fip05(:,i) = ode(xip05(i),yip05(:,i),FcnArgs{:});
        end
    end
    Fmid(:,iidx) = Fip05;
    F(:,xidx) = Freg;
end