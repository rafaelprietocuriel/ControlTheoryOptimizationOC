function [F Fmid Xmid Ymid]=calcdata_bvp6c(x,y,freepar,modelpar,ode)

global OCMATCONT OCBVP

% multi-point BVP support
arcposition = find(diff(x) == 0);
Lidx = [1, arcposition+1];
Ridx = [arcposition, length(x)];
FcnArgs = {0,freepar,modelpar};   
N=numel(x);

F = zeros(OCBVP.numode,N);
Fmid = zeros(OCBVP.numode,N-1);    % include interface intervals
for region = 1:OCBVP.numarc
    FcnArgs{1} = region;
    xidx = Lidx(region):Ridx(region);   % mesh point index
    Nreg = length(xidx);
    xreg = x(xidx);
    yreg = y(:,xidx);

    iidx = xidx(1:end-1);   % mesh interval index
    Nint = length(iidx);

    % derivative at the mesh points
    if OCMATCONT.OPTIONS.xyvectorized
        Freg = ode(xreg,yreg,FcnArgs{:});
    else
        Freg = zeros(OCBVP.numode,Nreg);
        for i = 1:Nreg
            Freg(:,i) = ode(xreg(i),yreg(:,i),FcnArgs{:});
        end
    end

    %mesh point data
    h = diff(xreg);
    H = spdiags(h(:),0,Nint,Nint);
    xi = xreg(1:end-1);
    yi = yreg(:,1:end-1);
    xip1 = xreg(2:end);
    yip1 = yreg(:,2:end);
    Fi = Freg(:,1:end-1);
    Fip1 = Freg(:,2:end);

    %interior points & derivative
    xip025 = 0.25*(3*xi + xip1);
    xip075 = 0.25*(xi + 3*xip1);
    yip025 = (54*yi + 10*yip1 + (9*Fi - 3*Fip1)*H)/64;
    yip075 = (10*yi + 54*yip1 + (3*Fi - 9*Fip1)*H)/64;
    if OCMATCONT.OPTIONS.xyvectorized
        Fip025 = ode(xip025,yip025,FcnArgs{:});
        Fip075 = ode(xip075,yip075,FcnArgs{:});
    else % not vectorized
        Fip025 = zeros(OCBVP.numode,Nint);
        Fip075 = zeros(OCBVP.numode,Nint);
        for i = 1:Nint
            Fip025(:,i) = ode(xip025(i),yip025(:,i),FcnArgs{:});
            Fip075(:,i) = ode(xip075(i),yip075(:,i),FcnArgs{:});
        end
    end

    %mid points & derivative
    xip05 = 0.5*(xi + xip1);
    yip05 = 0.5*(yi + yip1) - (5*Fi - 16*Fip025 + 16*Fip075 - 5*Fip1)*(H/24);
    if OCMATCONT.OPTIONS.xyvectorized
        Fip05 = ode(xip05,yip05,FcnArgs{:});
    else % not vectorized
        Fip05 = zeros(OCBVP.numode,Nint);
        for i = 1:Nint
            Fip05(:,i) = ode(xip05(i),yip05(:,i),FcnArgs{:});
        end
    end

    F(:,xidx) = Freg;
    Fmid(:,iidx,1) = Fip025;
    Fmid(:,iidx,2) = Fip05;
    Fmid(:,iidx,3) = Fip075;

    Xmid(:,iidx,1) = xip025;
    Xmid(:,iidx,2) = xip05;
    Xmid(:,iidx,3) = xip075;

    Ymid(:,iidx,1) = yip025;
    Ymid(:,iidx,2) = yip05;
    Ymid(:,iidx,3) = yip075;
end