function res=residual_bvp6c(x,y,freepar,modelpar,RHS,ode)
global OCMATCONT OCBVP

Fmid=OCBVP.Fmid;
yp=OCBVP.F;
FcnArgs = {0,freepar,modelpar};   

threshold=OCBVP.threshold;

lob(2)= 0.0848880518607165;
lobw(2)=0.276826047361566;
lob(3)= 0.265575603264643;
lobw(3)=0.431745381209863;

lob(5)= 0.734424396735357;
lobw(5)=lobw(3);
lob(6)= 0.9151119481392835;
lobw(6)=lobw(2);

Yp05=zeros(OCBVP.numode,OCBVP.N-1);    %output more accurate yp_ip05
Y05=zeros(OCBVP.numode,OCBVP.N-1);    %output more accurate yp_ip05

res = [];
nfcn = 0;

for arc = 1:OCBVP.numarc
    FcnArgs{1} = arc;    % Pass the arc index to the ODE function.

    xidx = OCBVP.Lidx(arc):OCBVP.Ridx(arc);    % mesh point index
    Nreg = OCBVP.Nint(arc)+1;
    xreg = x(xidx);
    yreg = y(:,xidx);
    ypreg = yp(:,xidx);
    hreg = diff(xreg);

    iidx = xidx(1:end-1);   % mesh interval index
    thresh = threshold(:,ones(1,OCBVP.Nint(arc)));
    

    yp_ip025=Fmid(:,iidx,1);
    yp_ip075=Fmid(:,iidx,3);

    diagscal = spdiags(hreg',0,OCBVP.Nint(arc),OCBVP.Nint(arc));
    xip05 = 0.5*(x(iidx) + x(iidx+1));
    %more accurate estimate of y_ip05 than used in Cash-Singhal
    Yip05 = 0.5*(y(:,iidx+1)+ y(:,iidx)) - ...
        (yp(:,iidx+1)-yp(:,iidx)+4*(yp_ip075-yp_ip025))*diagscal/24;
    if OCMATCONT.OPTIONS.xyvectorized
        Yp_ip05 = ode(xip05,Yip05,FcnArgs{:});
        nfcn = nfcn + 1;
    else % not vectorized
        Yp_ip05 = zeros(OCBVP.numode,OCBVP.Nint(arc));
        for i = 1:OCBVP.Nint(arc)
            Yp_ip05(:,i) = ode(xip05(i),Yip05(:,i),FcnArgs{:});
        end
        nfcn = nfcn + OCBVP.Nint(arc);
    end
    Y05(:,iidx)=Yip05;
    Yp05(:,iidx)=Yp_ip05;

    res_reg=zeros(1,OCBVP.Nint(arc));

    %Sum contributions from other points
    for j=[2,3,5,6]
        xLob = xreg(1:Nreg-1) + lob(j)*hreg;
        [yLob,ypLob] = interp_Hermite_bvp6c(lob(j),hreg,yreg,ypreg,yp_ip025,Yp_ip05,yp_ip075);
        if OCMATCONT.OPTIONS.xyvectorized
            fLob = ode(xLob,yLob,FcnArgs{:});
            nfcn = nfcn + 1;
        else
            fLob = zeros(OCBVP.numode,OCBVP.Nint(arc));
            for i = 1:OCBVP.Nint(arc)
                fLob(:,i) = ode(xLob(i),yLob(:,i),FcnArgs{:});
            end
            nfcn = nfcn + OCBVP.Nint(arc);
        end
        temp = (ypLob - fLob) ./ max(abs(fLob),thresh);
        res_reg = res_reg + lobw(j)*dot(temp,temp,1);
    end

    % scaling
    res_reg = sqrt( abs(hreg/2) .* res_reg);

    res(iidx) = res_reg;
end

OCBVP.Fmid(:,:,2)=Yp05;
OCBVP.Y05=Y05;

