function Jac=calcresjacpde_bvp4c(t,y,freepar,modelpar,ode,bc,odejac,bcjac)
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

%nfcn = 0;
%nbc = 0;

% BC points
ya = y(:,OCBVP.Lidx);
yb = y(:,OCBVP.Ridx);

rows = OCBVP.rows+OCBVP.numic;   % define the action area
cols = OCBVP.cols;             % in the global Jacobian
Jac = spalloc(OCBVP.nN+OCBVP.nparmc,OCBVP.nN+OCBVP.npar,(OCBVP.nN+OCBVP.npar)*(2*OCBVP.numode+OCBVP.npar));  % sparse storage
last_cols = zeros(OCBVP.nN+OCBVP.nparmc,OCBVP.npar);   % accumulator

if isempty(bcjac)   % use numerical approx
    [dGdya,dGdyb,nbc,dGdpar] = BCnumjac(bc,ya,yb,FcnArgs(2:3));
else  % use analytical Jacobian
    [dGdya,dGdyb,dGdpar] = bcjac(ya,yb,FcnArgs{2:3});
%     [dGdya1,dGdyb1,nbc,dGdpar1] = BCnumjac(bc,ya,yb,FcnArgs(2:3));
%     if any(abs(dGdya(:)-dGdya1(:))>1e-4) || any(abs(dGdyb(:)-dGdyb1(:))>1e-4) || any(abs(dGdpar(:)-dGdpar1(:))>1e-4)
%         dGdpar
%     end
end
last_cols(1:OCBVP.nBCs,:) = dGdpar;

% Collocation equations

for arc = 1:OCBVP.numarc

    % Left BC
    Jac(1:OCBVP.nBCs,cols) = dGdya(:,(arc-1)*OCBVP.numode+(1:OCBVP.numode));

    FcnArgs{1} = arc;

    tidx = OCBVP.Lidx(arc):OCBVP.Ridx(arc);
    treg = t(tidx);
    yreg = y(:,tidx);
    Freg = F(:,tidx);
    hreg = diff(treg);
    
    %iidx = tidx(1:end-1);    % mesh interval index
    Fmidreg = Fmid(:,tidx(1:end-1));
    % Collocation equations
    if isempty(odejac)  % use numerical approx
        OCBVP.Fref=Freg(:,1);
        [Ji,nFcalls,dFdpar_i]=Fnumjac(ode,{treg(1),yreg(:,1),FcnArgs{:}});
        %nfcn = nfcn+nFcalls;
        nrmJi = norm(Ji,1);
        nrmdFdpar_i = norm(dFdpar_i,1);

        for i = 1:OCBVP.Nint(arc)
            % the left mesh point
            ti = treg(i);
            yi = yreg(:,i);
            Fi = Freg(:,i);
            % the right mesh point
            tip1 = treg(i+1);
            yip1 = yreg(:,i+1);
            Fip1 = Freg(:,i+1);
            OCBVP.Fref=Fip1;
            [Jip1,nFcalls,dFdpar_ip1]=Fnumjac(ode,{tip1,yip1,FcnArgs{:}});
            %nfcn = nfcn + nFcalls;
            nrmJip1 = norm(Jip1,1);
            nrmdFdpar_ip1 = norm(dFdpar_ip1,1);
            % the midpoint
            hi = hreg(i);
            tip05 = (ti + tip1)/2;
            if (norm(Jip1 - Ji,1) <= 0.25*(nrmJi + nrmJip1)) && ...
                    (norm(dFdpar_ip1 - dFdpar_i,1) <= 0.25*(nrmdFdpar_i + nrmdFdpar_ip1))
                Jip05 = 0.5*(Ji + Jip1);
                dFdpar_ip05 = 0.5*(dFdpar_i + dFdpar_ip1);
            else
                yip05 = (yi+yip1)/2-hi/8*(Fip1-Fi);
                OCBVP.Fref=Fmidreg(:,i);
                [Jip05,nFcalls,dFdpar_ip05]=Fnumjac(ode,{tip05,yip05,FcnArgs{:}});
                %nfcn = nfcn + nFcalls;1
            end
            twiceJip05 = 2*Jip05;
            % assembly
            Jac(rows,cols) = -(OCBVP.In+hi/6*(Ji+twiceJip05*(OCBVP.In+hi/4*Ji)));
            cols = cols + OCBVP.numode;
            Jac(rows,cols) = OCBVP.In-hi/6*(Jip1+twiceJip05*(OCBVP.In-hi/4*Jip1));
            last_cols(rows,:) = -hi*dFdpar_ip05 + hi^2/12*Jip05*...
                (dFdpar_ip1-dFdpar_i);
            rows = rows+OCBVP.numode;   % next equation

            Ji = Jip1;
            nrmJi = nrmJip1;
            dFdpar_i = dFdpar_ip1;
            nrmdFdpar_i = nrmdFdpar_ip1;
        end

    else % use analytical Jacobian

        [Ji,dFdpar_i] = odejac(treg(1),yreg(:,1),FcnArgs{:});
%         OCBVP.Fref=Freg(:,1);
%         [Jitest,nFcalls,dFdpar_itest]=Fnumjac(ode,{treg(1),yreg(:,1),FcnArgs{:}});
%         if any(abs(Ji(:)-Jitest(:))>1e-1) || any(abs(dFdpar_i(:)-dFdpar_itest(:))>1e0)
%             max(abs(dFdpar_i(:)-dFdpar_itest(:)))
%             max(abs(Ji(:)-Jitest(:)))
%         end


        for i = 1:OCBVP.Nint(arc)
            % the left mesh point
            ti = treg(i);
            yi = yreg(:,i);
            Fi = Freg(:,i);
            % the right mesh point
            tip1 = treg(i+1);
            yip1 = yreg(:,i+1);
            Fip1 = Freg(:,i+1);
            [Jip1, dFdpar_ip1] = odejac(tip1,yip1,FcnArgs{:});
            % the midpoint
            hi = hreg(i);
            tip05 = (ti + tip1)/2;
            yip05 = (yi + yip1)/2 - hi/8*(Fip1 - Fi);
            [Jip05, dFdpar_ip05] = odejac(tip05,yip05,FcnArgs{:});  % recompute the Jacobian
            % test
            %yip05 = (yi+yip1)/2-hi/8*(Fip1-Fi);
%             OCBVP.Fref=Fmidreg(:,i);
%             [Jip05test,nFcalls,dFdpar_ip05test]=Fnumjac(ode,{tip05,yip05,FcnArgs{:}});
%             if any(abs(Jip05(:)-Jip05test(:))>1e-1) || any(abs(dFdpar_ip05(:)-dFdpar_ip05test(:))>1e-1)
%                 max(abs(Jip05(:)-Jip05test(:)))
%             end

            twiceJip05 = 2*Jip05;
            % assembly
            Jac(rows,cols) = -(OCBVP.In+hi/6*(Ji+twiceJip05*(OCBVP.In+hi/4*Ji)));
            cols = cols + OCBVP.numode;
            Jac(rows,cols) = OCBVP.In-hi/6*(Jip1+twiceJip05*(OCBVP.In-hi/4*Jip1));
            last_cols(rows,:) = -hi*dFdpar_ip05 + hi^2/12*Jip05* ...
                (dFdpar_ip1-dFdpar_i);
            rows = rows+OCBVP.numode;   % next equation

            Ji = Jip1;
            dFdpar_i = dFdpar_ip1;
        end
    end

    % Right BC
    Jac(1:OCBVP.nBCs,cols) = dGdyb(:,(arc-1)*OCBVP.numode+(1:OCBVP.numode));
    cols = cols + OCBVP.numode;
end

Jac(:,end-OCBVP.npar+1:end) = last_cols;  % accumulated

