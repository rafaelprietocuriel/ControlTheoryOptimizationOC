function M=calcmonodromy_bvp6c(x,y,freepar,modelpar,ode,odejac,F,Fmid,Xmid,Ymid)

global OCBVP
FcnArgs = {0,freepar,modelpar};    % Pass the arc index to the ODE function.

threshval = 1e-6;
OCBVP.Joptions.thresh = threshval(ones(OCBVP.numode,1));
OCBVP.dPoptions.thresh = threshval(ones(OCBVP.npar,1));

M=eye(OCBVP.numode);

averageJac=OCBVP.averageJac;

rows = OCBVP.rows;   % define the action area
cols = OCBVP.cols;             % in the global Jacobian
counter=0;
% Collocation equations
for region = 1:OCBVP.numarc

    FcnArgs{1} = region;

    xidx = OCBVP.Lidx(region):OCBVP.Ridx(region);
    xreg = x(xidx);
    yreg = y(:,xidx);
    hreg = diff(xreg);

    iidx = xidx(1:end-1);    % mesh interval index
    Nint = length(iidx);
        Freg = F(:,xidx);
        hreg = diff(xreg);

    [X1qtrreg, Xmidreg, X3qtrreg] = midptreg(iidx,Xmid);
    [Y1qtrreg, Ymidreg, Y3qtrreg] = midptreg(iidx,Ymid);
    [F1qtrreg, Fmidreg] = midptreg(iidx,Fmid);
    % Collocation equations
    if isempty(odejac)  % use numerical approx
        OCBVP.Fref=Freg(:,1);
        [Ji,nFcalls]=Fnumjac(ode,{xreg(1),yreg(:,1),FcnArgs{:}});
        nfcn = nfcn+nFcalls;
        nrmJi = norm(Ji,1);

        for i = 1:Nint
            counter=counter+1;
            hi = hreg(i);
            % the right mesh point
            xip1 = xreg(i+1);
            yip1 = yreg(:,i+1);
            Fip1 = Freg(:,i+1);
            OCBVP.Fref=Fip1;
            [Jip1]=Fnumjac(ode,{xip1,yip1,FcnArgs{:}});
            nrmJip1 = norm(Jip1,1);

            %the interior points
            if averageJac && (norm(Jip1 - Ji,1) <= 0.125*(nrmJi + nrmJip1))
                Jip025 = 0.25*(3*Ji + Jip1);
                Jip05 = 0.5*(Ji + Jip1);
                Jip075 = 0.25*(Ji + 3*Jip1);
            else
                [xip025, xip05, xip075] = midpti(i,X1qtrreg, Xmidreg, X3qtrreg);
                [yip025, yip05, yip075] = midpti(i,Y1qtrreg, Ymidreg, Y3qtrreg);

                OCBVP.Fref=Fmidreg(:,i);
                [Jip025]=Fnumjac(ode,{xip025,yip025,FcnArgs{:}});
                [Jip05]=Fnumjac(ode,{xip05,yip05,FcnArgs{:}});
                [Jip075]=Fnumjac(ode,{xip075,yip075,FcnArgs{:}});
            end

            Jip05Jip025=Jip05*Jip025;
            Jip05Jip075=Jip05*Jip075;
            % assembly
            Gi=calc_dPHYdy1(hi,OCBVP.In,Ji,Jip025,Jip05,Jip075,Jip05Jip025,Jip05Jip075);
            cols = cols + OCBVP.numode;
            Hi=calc_dPHYdy2(hi,OCBVP.In,Jip1,Jip025,Jip05,Jip075,Jip05Jip025,Jip05Jip075);
            M=-inv(Hi)*Gi*M;
            rows = rows + OCBVP.numode;   % next equation

            Ji = Jip1;
            nrmJi = nrmJip1;
            OCBVP.Hi(:,:,counter)=Hi;
            OCBVP.Gi(:,:,counter)=Gi;
        end

    else % use analytical Jacobian

        [Ji] = odejac(xreg(1),yreg(:,1),FcnArgs{:});

        for i = 1:Nint
            counter=counter+1;
            hi = hreg(i);
            % the right mesh point
            xip1 = xreg(i+1);
            yip1 = yreg(:,i+1);
            [Jip1] = odejac(xip1,yip1,FcnArgs{:});

            %the interior points
            [xip025, xip05, xip075] = midpti(i,X1qtrreg, Xmidreg, X3qtrreg);
            [yip025, yip05, yip075] = midpti(i,Y1qtrreg, Ymidreg, Y3qtrreg);

            [Jip025] = odejac(xip025,yip025,FcnArgs{:});
            [Jip05] = odejac(xip05,yip05,FcnArgs{:});
            [Jip075] = odejac(xip075,yip075,FcnArgs{:});

            Jip05Jip025=Jip05*Jip025;
            Jip05Jip075=Jip05*Jip075;
            % assembly
            Gi=calc_dPHYdy1(hi,OCBVP.In,Ji,Jip025,Jip05,Jip075,Jip05Jip025,Jip05Jip075);
            cols = cols + OCBVP.numode;
            Hi=calc_dPHYdy2(hi,OCBVP.In,Jip1,Jip025,Jip05,Jip075,Jip05Jip025,Jip05Jip075);
            %M=-inv(Hi)*Gi*M;
            M=-Hi\Gi*M;
            rows = rows + OCBVP.numode;   % next equation

            Ji = Jip1;
            OCBVP.Hi(:,:,counter)=Hi;
            OCBVP.Gi(:,:,counter)=Gi;
        end
    end

    cols = cols + OCBVP.numode;
end

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


