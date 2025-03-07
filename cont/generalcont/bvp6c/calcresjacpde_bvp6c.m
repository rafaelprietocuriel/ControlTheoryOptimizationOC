function Jac=calcresjacpde_bvp6c(t,y,freepar,modelpar,ode,bc,odejac,bcjac)

global OCBVP
FcnArgs = {0,freepar,modelpar};    % Pass the arc index to the ODE function.

threshval = 1e-6;
OCBVP.Joptions.thresh = threshval(ones(OCBVP.numode,1));
OCBVP.dPoptions.thresh = threshval(ones(OCBVP.npar,1));

F=OCBVP.F;
Fmid=OCBVP.Fmid;
Xmid=OCBVP.Xmid;
Ymid=OCBVP.Ymid;

averageJac=OCBVP.averageJac;
nfcn = 0;

% BC points
ya = y(:,OCBVP.Lidx);
yb = y(:,OCBVP.Ridx);

rows = OCBVP.rows+OCBVP.numic;   % define the action area
cols = OCBVP.cols;             % in the global Jacobian
icolsstart=0;

Jac = spalloc(OCBVP.nN+OCBVP.nparmc,OCBVP.nN+OCBVP.npar,(OCBVP.nN+OCBVP.npar)*(2*OCBVP.numode+OCBVP.npar));  % sparse storage
last_cols = zeros(OCBVP.nN+OCBVP.nparmc,OCBVP.npar);   % accumulator

if isempty(bcjac)   % use numerical approx
    [dGdya,dGdyb,nbc,dGdpar] = BCnumjac(bc,ya,yb,{freepar,modelpar});
else  % use analytical Jacobian
    [dGdya,dGdyb,dGdpar]=bcjac(ya,yb,freepar,modelpar);
end
last_cols(1:OCBVP.nBCs,:) = dGdpar;

% Collocation equations
for region = 1:OCBVP.numarc

    % Left BC
    Jac(1:OCBVP.nBCs,cols) = dGdya(:,(region-1)*OCBVP.numode+(1:OCBVP.numode));

    FcnArgs{1} = region;

    tidx = OCBVP.Lidx(region):OCBVP.Ridx(region);
    treg = t(tidx);
    yreg = y(:,tidx);
    Freg = F(:,tidx);
    hreg = diff(treg);

    iidx = tidx(1:end-1);    % mesh interval index
    Nint = length(iidx);

    [X1qtrreg, Xmidreg, X3qtrreg] = midptreg(iidx,Xmid);
    [Y1qtrreg, Ymidreg, Y3qtrreg] = midptreg(iidx,Ymid);
    [F1qtrreg, Fmidreg, F3qtrreg] = midptreg(iidx,Fmid);

    % Collocation equations
    if isempty(odejac)  % use numerical approx
        OCBVP.Fref=Freg(:,1);
        [Ji,nFcalls,dFdpar_i]=Fnumjac(ode,{treg(1),yreg(:,1),FcnArgs{:}});
        nfcn = nfcn+nFcalls;
        nrmJi = norm(Ji,1);
        nrmdFdpar_i = norm(dFdpar_i,1);

        for i = 1:Nint
            hi = hreg(i);
            % the right mesh point
            tip1 = treg(i+1);
            yip1 = yreg(:,i+1);
            Fip1 = Freg(:,i+1);
            OCBVP.Fref=Fip1;
            [Jip1,nFcalls,dFdpar_ip1]=Fnumjac(ode,{tip1,yip1,FcnArgs{:}});
            nfcn = nfcn + nFcalls;
            nrmJip1 = norm(Jip1,1);
            nrmdFdpar_ip1 = norm(dFdpar_ip1,1);

            %the interior points
            if averageJac && (norm(Jip1 - Ji,1) <= 0.125*(nrmJi + nrmJip1)) && ...
                    (norm(dFdpar_ip1 - dFdpar_i,1) <= 0.125*(nrmdFdpar_i + nrmdFdpar_ip1))
                Jip025 = 0.25*(3*Ji + Jip1);
                Jip05 = 0.5*(Ji + Jip1);
                Jip075 = 0.25*(Ji + 3*Jip1);

                dFdpar_ip025 = 0.25*(3*dFdpar_i + dFdpar_ip1);
                dFdpar_ip05 = 0.5*(dFdpar_i + dFdpar_ip1);
                dFdpar_ip075 = 0.25*(dFdpar_i + 3*dFdpar_ip1);
            else
                [tip025, tip05, tip075] = midpti(i,X1qtrreg, Xmidreg, X3qtrreg);
                [yip025, yip05, yip075] = midpti(i,Y1qtrreg, Ymidreg, Y3qtrreg);

                OCBVP.Fref=Fmidreg(:,i);
                [Jip025,nFcalls025,dFdpar_ip025]=Fnumjac(ode,{tip025,yip025,FcnArgs{:}});
                [Jip05,nFcalls05,dFdpar_ip05]=Fnumjac(ode,{tip05,yip05,FcnArgs{:}});
                [Jip075,nFcalls075,dFdpar_ip075]=Fnumjac(ode,{tip075,yip075,FcnArgs{:}});

                nfcn = nfcn + nFcalls025 + nFcalls05 + nFcalls075;
            end

            Jip05Jip025=Jip05*Jip025;
            Jip05Jip075=Jip05*Jip075;
            % assembly
            Jac(rows,cols)=calc_dPHYdy1(hi,OCBVP.In,Ji,Jip025,Jip05,Jip075,Jip05Jip025,Jip05Jip075);
            cols = cols + OCBVP.numode;
            Jac(rows,cols)=calc_dPHYdy2(hi,OCBVP.In,Jip1,Jip025,Jip05,Jip075,Jip05Jip025,Jip05Jip075);

            last_cols(rows,:) = - ( ...
                hi/90*( 7*(dFdpar_i+dFdpar_ip1)+32*(dFdpar_ip025+dFdpar_ip075)+12*dFdpar_ip05 )     + ...
                hi*hi/180*( Jip025*(9*dFdpar_i-3*dFdpar_ip1) + Jip075*(3*dFdpar_i-9*dFdpar_ip1) + ...
                Jip05*( 16*(dFdpar_ip025-dFdpar_ip075)+5*(dFdpar_ip1-dFdpar_i) )  )     + ...
                hi*hi*hi/240*( Jip05Jip025*(3*dFdpar_i-dFdpar_ip1)+Jip05Jip075*(3*dFdpar_ip1-dFdpar_i) )  ...
                );
            rows = rows + OCBVP.numode;   % next equation

            Ji = Jip1;
            nrmJi = nrmJip1;
            dFdpar_i = dFdpar_ip1;
            nrmdFdpar_i = nrmdFdpar_ip1;
        end

    else % use analytical Jacobian

        [Ji,dFdpar_i] = odejac(treg(1),yreg(:,1),FcnArgs{:});

        for i = 1:Nint
            hi = hreg(i);
            % the right mesh point
            tip1 = treg(i+1);
            yip1 = yreg(:,i+1);
            [Jip1, dFdpar_ip1] = odejac(tip1,yip1,FcnArgs{:});

            %the interior points
            [tip025, tip05, tip075] = midpti(i,X1qtrreg, Xmidreg, X3qtrreg);
            [yip025, yip05, yip075] = midpti(i,Y1qtrreg, Ymidreg, Y3qtrreg);

            [Jip025, dFdpar_ip025] = odejac(tip025,yip025,FcnArgs{:});
            [Jip05, dFdpar_ip05] = odejac(tip05,yip05,FcnArgs{:});
            [Jip075, dFdpar_ip075] = odejac(tip075,yip075,FcnArgs{:});

            Jip05Jip025=Jip05*Jip025;
            Jip05Jip075=Jip05*Jip075;
            % assembly
            Jac(rows,cols)=calc_dPHYdy1(hi,OCBVP.In,Ji,Jip025,Jip05,Jip075,Jip05Jip025,Jip05Jip075);
            cols = cols + OCBVP.numode;
            Jac(rows,cols)=calc_dPHYdy2(hi,OCBVP.In,Jip1,Jip025,Jip05,Jip075,Jip05Jip025,Jip05Jip075);

            last_cols(rows,:) = - ( ...
                hi/90*( 7*(dFdpar_i+dFdpar_ip1)+32*(dFdpar_ip025+dFdpar_ip075)+12*dFdpar_ip05 )     + ...
                hi*hi/180*( Jip025*(9*dFdpar_i-3*dFdpar_ip1) + Jip075*(3*dFdpar_i-9*dFdpar_ip1) + ...
                Jip05*( 16*(dFdpar_ip025-dFdpar_ip075)+5*(dFdpar_ip1-dFdpar_i) )  )     + ...
                hi*hi*hi/240*( Jip05Jip025*(3*dFdpar_i-dFdpar_ip1)+Jip05Jip075*(3*dFdpar_ip1-dFdpar_i) )  ...
                );
            rows = rows + OCBVP.numode;   % next equation

            Ji = Jip1;
            dFdpar_i = dFdpar_ip1;
        end
    end

    % Right BC
    Jac(1:OCBVP.nBCs,cols) = dGdyb(:,(region-1)*OCBVP.numode+(1:OCBVP.numode));

    cols = cols + OCBVP.numode;
end

Jac(:,end-OCBVP.npar+1:end) = last_cols;  % accumulated

%--------------------------------------------------------------------------

function J=calc_dPHYdy1(hi,In,Ji,Jip025,Jip05,Jip075,Jip05Jip025,Jip05Jip075)

J   = - In -             ( ...
    hi/90*       ( 7*Ji+27*Jip025+6*Jip05+5*Jip075   ) + ...
    hi*hi/360*   ( 27*Jip05Jip025-5*Jip05Jip075)       + ...
    ( hi*hi/360*(18*Jip025-10*Jip05+6*Jip075) + ...
    hi*hi*hi/240*(3*Jip05Jip025-Jip05Jip075) )*Ji      ...
    );
%--------------------------------------------------------------------------

function J=calc_dPHYdy2(hi,In,Jip1,Jip025,Jip05,Jip075,Jip05Jip025,Jip05Jip075)

J   = In -               (...
    hi/90*       ( 5*Jip025+6*Jip05+27*Jip075+7*Jip1  ) + ...
    hi*hi/360*   ( 5*Jip05Jip025-27*Jip05Jip075 )       - ...
    ( hi*hi/360*(6*Jip025-10*Jip05+18*Jip075) + ...
    hi*hi*hi/240*(Jip05Jip025-3*Jip05Jip075)  )*Jip1   ...
    );


%--------------------------------------------------------------------------

function [a025, a05, a075] = midptreg(iidx,mid)
a025 = mid(:,iidx,1);
a05 = mid(:,iidx,2);
a075 = mid(:,iidx,3);

function [a025, a05, a075] = midpti(i, b025, b05, b075)
a025 = b025(:,i);
a05 = b05(:,i);
a075 = b075(:,i);


