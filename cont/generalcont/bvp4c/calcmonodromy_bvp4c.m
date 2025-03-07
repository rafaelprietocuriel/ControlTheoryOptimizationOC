function M=calcmonodromy_bvp4c(x,y,freepar,modelpar,ode,odejac,F,Fmid)
% This is a straight forward implementation for the computation of the
% monodromy matrix. Its simplicitiy lacks exactness if the eigenvalues of
% the Floquet multiplier differ substantially in size. To overcome this
% problem the periodic Schur decomposition could be applied as described in
% lust1997, lust2001
%
% @PHDTHESIS{lust1997,
%   author = {Lust, K.},
%   title = {Numerical Bifurcation Analysis of Periodic Solutions of Partial Differential
% 	Equations},
%   school = {Department of Computer Science, K.U.Leuven},
%   year = {1997}}
%
% @ARTICLE{lust2001,
%   author = {Lust, K.},
%   title = {Improved numerical Floquet multipliers},
%   journal = {International Journal of Bifurcation and Chaos},
%   year = {2001},
%   volume = {11},
%   pages = {2389--2410},
%   number = {9}}

global OCMATCONT OCBVP
FcnArgs = {0,freepar,modelpar};    % Pass the region index to the ODE function.

threshval = 1e-6;
OCBVP.Joptions.thresh = threshval(ones(OCBVP.numode,1));
OCBVP.dPoptions.thresh = threshval(ones(OCBVP.npar,1));

M=eye(OCBVP.numode);
rows = OCBVP.rows;   % define the action area
cols = OCBVP.cols;             % in the global Jacobian

OCBVP.Hi=zeros(OCBVP.numode,OCBVP.numode,length(x));
OCBVP.Gi=zeros(OCBVP.numode,OCBVP.numode,length(x));

counter=0;
% Collocation equations
for arc = 1:OCBVP.numarc

    FcnArgs{1} = arc;

    xidx = OCBVP.Lidx(arc):OCBVP.Ridx(arc);
    xreg = x(xidx);
    yreg = y(:,xidx);
    hreg = diff(xreg);
    xi = xreg(1:end-1);
    yi = yreg(:,1:end-1);
    xip1 = xreg(2:end);
    yip1 = yreg(:,2:end);
    Freg = F(:,xidx);
    Fmidreg = Fmid(:,xidx(1:end-1));
    %iidx = xidx(1:end-1);    % mesh interval index

    % Collocation equations
    if isempty(odejac)  % use numerical approx
        OCBVP.Fref=Freg(:,1);
        [Ji,nFcalls,dFdpar_i]=Fnumjac(ode,{xreg(1),yreg(:,1),FcnArgs{:}});
        %nfcn = nfcn+nFcalls;
        nrmJi = norm(Ji,1);
        nrmdFdpar_i = norm(dFdpar_i,1);

        for i = 1:OCBVP.Nint(arc)
            counter=counter+1;
            % the left mesh point
            xi = xreg(i);
            yi = yreg(:,i);
            Fi = Freg(:,i);
            % the right mesh point
            xip1 = xreg(i+1);
            yip1 = yreg(:,i+1);
            Fip1 = Freg(:,i+1);
            OCBVP.Fref=Fip1;
            [Jip1,nFcalls,dFdpar_ip1]=Fnumjac(ode,{xip1,yip1,FcnArgs{:}});
            %nfcn = nfcn + nFcalls;
            nrmJip1 = norm(Jip1,1);
            nrmdFdpar_ip1 = norm(dFdpar_ip1,1);
            % the midpoint
            hi = hreg(i);
            xip05 = (xi + xip1)/2;
            if (norm(Jip1 - Ji,1) <= 0.25*(nrmJi + nrmJip1)) && ...
                    (norm(dFdpar_ip1 - dFdpar_i,1) <= 0.25*(nrmdFdpar_i + nrmdFdpar_ip1))
                Jip05 = 0.5*(Ji + Jip1);
            else
                yip05 = (yi+yip1)/2-hi/8*(Fip1-Fi);
                OCBVP.Fref=Fmidreg(:,i);
                Jip05=Fnumjac(ode,{xip05,yip05,FcnArgs{:}});
                %nfcn = nfcn + nFcalls;
            end
            twiceJip05 = 2*Jip05;
            % assembly
            Gi=-(OCBVP.In+hi/6*(Ji+twiceJip05*(OCBVP.In+hi/4*Ji)));
            cols = cols + OCBVP.numode;
            Hi=OCBVP.In-hi/6*(Jip1+twiceJip05*(OCBVP.In-hi/4*Jip1));
            M=-inv(Hi)*Gi*M;
            rows = rows+OCBVP.numode;   % next equation

            Ji = Jip1;
            nrmJi = nrmJip1;
            dFdpar_i = dFdpar_ip1;
            nrmdFdpar_i = nrmdFdpar_ip1;
            OCBVP.Hi(:,:,counter)=Hi;
            OCBVP.Gi(:,:,counter)=Gi;
        end

    else % use analytical Jacobian

        Ji=odejac(xreg(1),yreg(:,1),FcnArgs{:});

        for i = 1:OCBVP.Nint(arc)
            counter=counter+1;
            % the left mesh point
            xi = xreg(i);
            yi = yreg(:,i);
            Fi = Freg(:,i);
            % the right mesh point
            xip1 = xreg(i+1);
            yip1 = yreg(:,i+1);
            Fip1 = Freg(:,i+1);
            Jip1=odejac(xip1,yip1,FcnArgs{:});
            % the midpoint
            hi = hreg(i);
            xip05 = (xi + xip1)/2;
            yip05 = (yi + yip1)/2 - hi/8*(Fip1 - Fi);
            Jip05=odejac(xip05,yip05,FcnArgs{:});  % recompute the Jacobian
            twiceJip05 = 2*Jip05;
            % assembly
            Gi=-(OCBVP.In+hi/6*(Ji+twiceJip05*(OCBVP.In+hi/4*Ji)));
            cols = cols + OCBVP.numode;
            Hi=OCBVP.In-hi/6*(Jip1+twiceJip05*(OCBVP.In-hi/4*Jip1));
            M=-inv(Hi)*Gi*M;
            rows = rows+OCBVP.numode;   % next equation

            Ji = Jip1;
            OCBVP.Hi(:,:,counter)=Hi;
            OCBVP.Gi(:,:,counter)=Gi;
        end
    end
    cols = cols + OCBVP.numode;
end

