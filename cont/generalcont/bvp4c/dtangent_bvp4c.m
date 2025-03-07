function [Dv,Dvmid]=dtangent_bvp4c(x,y,v,freepar,modelpar,ode,odejac)
%

% This function bases on bvp4c by
%   Jacek Kierzenka and Lawrence F. Shampine
%   Copyright 1984-2007 The MathWorks, Inc.
%   $Revision: 1.21.4.15 $  $Date: 2007/05/23 18:54:07 $

global OCBVP
FcnArgs = {0,freepar,modelpar};    % Pass the region index to the ODE function.
threshval = 1e-6;
OCBVP.Joptions.thresh = threshval(ones(OCBVP.numode,1));
OCBVP.dPoptions.thresh = threshval(ones(OCBVP.npar,1));

F=OCBVP.F;
Fmid=OCBVP.Fmid;
Dv=zeros(size(F));
Dvmid=zeros(size(Fmid));
cols=1;
for arc = 1:OCBVP.numarc

    % Left BC
    FcnArgs{1} = arc;

    xidx = OCBVP.Lidxold(arc):OCBVP.Ridxold(arc);
    xreg = x(xidx);
    yreg = y(:,xidx);
    Freg = F(:,xidx);
    hreg = diff(xreg);
    vreg=v(:,xidx);
    %iidx = xidx(1:end-1);    % mesh interval index
    Fmidreg = Fmid(:,xidx(1:end-1));
    % Collocation equations
    if isempty(odejac)  % use numerical approx
        OCBVP.Fref=Freg(:,1);
        Ji=Fnumjac(ode,{xreg(1),yreg(:,1),FcnArgs{:}});
        Dv(:,cols)=Ji*vreg(:,1);

        for i = 1:OCBVP.Nintold(arc)
            % the left mesh point
            xi = xreg(i);
            yi = yreg(:,i);
            Fi = Freg(:,i);
            vi = vreg(:,i);
            % the right mesh point
            xip1 = xreg(i+1);
            yip1 = yreg(:,i+1);
            Fip1 = Freg(:,i+1);
            vip1 = vreg(:,i+1);
            OCBVP.Fref=Fip1;
            Jip1=Fnumjac(ode,{xip1,yip1,FcnArgs{:}});

            % the midpoint
            hi = hreg(i);
            xip05 = (xi + xip1)/2;
            yip05 = (yi+yip1)/2-hi/8*(Fip1-Fi);
            vip05 = (vi + vip1)/2;
            OCBVP.Fref=Fmidreg(:,i);
            Jip05=Fnumjac(ode,{xip05,yip05,FcnArgs{:}});

            % assembly
            Dvmid(:,cols) = Jip05*vip05;
            cols = cols + 1;
            Dv(:,cols) =Jip1*vip1;
        end
    else % use analytical Jacobian

        Ji=odejac(xreg(1),yreg(:,1),FcnArgs{:});
        Dv(:,cols)=Ji*vreg(:,1);

        for i = 1:OCBVP.Nintold(arc)
            % the left mesh point
            xi = xreg(i);
            yi = yreg(:,i);
            Fi = Freg(:,i);
            vi = vreg(:,i);
            % the right mesh point
            xip1 = xreg(i+1);
            yip1 = yreg(:,i+1);
            Fip1 = Freg(:,i+1);
            vip1 = vreg(:,i+1);
            Jip1=odejac(xip1,yip1,FcnArgs{:});
            % the midpoint
            hi = hreg(i);
            xip05 = (xi + xip1)/2;
            yip05 = (yi + yip1)/2 - hi/8*(Fip1 - Fi);
            vip05 = (vi + vip1)/2;
            Jip05=odejac(xip05,yip05,FcnArgs{:});  % recompute the Jacobian
            % assembly
            Dvmid(:,cols) = Jip05*vip05;
            cols = cols + 1;
            Dv(:,cols) =Jip1*vip1;
        end
    end
end
