    function Jac=calcresjac_gbvp4c(x,y,freepar,modelpar,ode,bc,icfun,odejac,bcjac,icfunjac)
%

% This function bases on bvp4c by
%   Jacek Kierzenka and Lawrence F. Shampine
%   Copyright 1984-2007 The MathWorks, Inc.
%   $Revision: 1.21.4.15 $  $Date: 2007/05/23 18:54:07 $

global OCMATCONT OCBVP
FcnArgs = {0,freepar,modelpar};    % Pass the region index to the ODE function.

F=OCBVP.F;
Fmid=OCBVP.Fmid;

%nfcn = 0;
%nbc = 0;

% BC points
ya = y(:,OCBVP.Lidx);
yb = y(:,OCBVP.Ridx);

rows = OCBVP.rows;   % define the action area
cols = 1:OCBVP.numode(1);             % in the global Jacobian

Jac = spalloc(OCBVP.nN+OCBVP.nparmc,OCBVP.nN+OCBVP.npar,(OCBVP.nN+OCBVP.npar)*(2*OCBVP.maxnumode+OCBVP.npar));  % sparse storage
last_cols = zeros(OCBVP.nN+OCBVP.nparmc,OCBVP.npar);   % accumulator

if OCMATCONT.OPTIONS.messagedisplay
    StartTime = clock;
end
if isempty(bcjac)   % use numerical approx
    [dGdya,dGdyb,nbc,dGdpar] = gBCnumjac(bc,ya,yb,FcnArgs(2:3));
    %     numJacOpt.diffvar=1;
    %     numJacOpt.vectvars=[];
    %     dGdya1test=numjaccsd4bc(bc,{ya(:),yb,FcnArgs{2:3}},OCBVP.nBCs,size(ya),numJacOpt);
    %     numJacOpt.diffvar=2;
    %     dGdyb1test=numjaccsd4bc(bc,{ya,yb(:),FcnArgs{2:3}},OCBVP.nBCs,size(ya),numJacOpt);
    %     numJacOpt.diffvar=3;
    %     dGdpartest=numjaccsd4bc(bc,{ya,yb,FcnArgs{2:3}},OCBVP.nBCs,size(freepar),numJacOpt);
else  % use analytical Jacobian
    [dGdya,dGdyb,dGdpar] = bcjac(ya,yb,FcnArgs{2:3});
    %[dGdya1,dGdyb1,nbc,dGdpar1] = BCnumjac(bc,ya,yb,FcnArgs(2:3));
    %     if any(abs(dGdya(:)-dGdya1(:))>1e-4) || any(abs(dGdyb(:)-dGdyb1(:))>1e-4) || any(abs(dGdpar(:)-dGdpar1(:))>1e-4)
    %         dGdpar
    %     end
end
last_cols(1:OCBVP.nBCs,:) = dGdpar;
if OCMATCONT.OPTIONS.messagedisplay
    ocmatmsg('elapsed time (DG)=%.1f secs\n', etime(clock, StartTime));
    StartTime = clock;
end

for arc = 1:OCBVP.numarc
    OCBVP.neqn=OCBVP.numode(arc);
    % Left BC
    Jac(1:OCBVP.nBCs,cols) = dGdya(:,sum(OCBVP.numode(1:arc-1))+(1:OCBVP.numode(arc)));

    FcnArgs{1} = arc;

    xidx = OCBVP.Lidx(arc):OCBVP.Ridx(arc);
    xreg = x(xidx);
    yreg = y(1:OCBVP.numode(arc),xidx);
    Freg = F(1:OCBVP.numode(arc),xidx);
    hreg = diff(xreg);
    %iidx = xidx(1:end-1);    % mesh interval index
    Fmidreg = Fmid(:,xidx(1:end-1));

    % Collocation equations
    if isempty(odejac)  % use numerical approx
        OCBVP.Fref=Freg(:,1);
        [Ji,nFcalls,dFdpar_i]=Fnumjac(ode,{xreg(1),yreg(:,1),FcnArgs{:}});

        nrmJi = norm(Ji,1);
        nrmdFdpar_i = norm(dFdpar_i,1);

        for i = 1:OCBVP.Nint(arc)
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

            nrmJip1 = norm(Jip1,1);
            nrmdFdpar_ip1 = norm(dFdpar_ip1,1);
            % the midpoint
            hi = hreg(i);
            xip05 = (xi + xip1)/2;
            if (norm(Jip1 - Ji,1) <= 0.25*(nrmJi + nrmJip1)) && ...
                    (norm(dFdpar_ip1 - dFdpar_i,1) <= 0.25*(nrmdFdpar_i + nrmdFdpar_ip1))
                Jip05 = 0.5*(Ji + Jip1);
                dFdpar_ip05 = 0.5*(dFdpar_i + dFdpar_ip1);
            else
                yip05 = (yi+yip1)/2-hi/8*(Fip1-Fi);
                OCBVP.Fref=Fmidreg(:,i);
                [Jip05,nFcalls,dFdpar_ip05]=Fnumjac(ode,{xip05,yip05,FcnArgs{:}});
            end

            twiceJip05 = 2*Jip05;
            % assembly
            Jac(rows,cols) = -(OCBVP.In{arc}+hi/6*(Ji+twiceJip05*(OCBVP.In{arc}+hi/4*Ji)));
            cols = cols + OCBVP.numode(arc);
            Jac(rows,cols) = OCBVP.In{arc}-hi/6*(Jip1+twiceJip05*(OCBVP.In{arc}-hi/4*Jip1));
            % the exact local Jacobian with respect to par is
            % -1/6*(dFdpar_ip1-1/2*Jip05*(dFdpar_ip1-dFdpar_i)*hi+4*dFdpar_ip05+dFdpar_i)*hi
            % for the actual computation the approximation (dFdpar_ip1+dFdpar_i)/2~dFdpar_ip05 is used.
            % therefore: dFdpar_ip1+dFdpar_i+4*dFdpar_ip05~6*dFdpar_ip05
            % simplifies to

            %             last_cols(rows,:) = -hi*dFdpar_ip05 + hi^2/12*Jip05*...
            %                 (dFdpar_ip1-dFdpar_i);
            last_cols(rows,:) = -1/6*(dFdpar_ip1-1/2*Jip05*(dFdpar_ip1-dFdpar_i)*hi+4*dFdpar_ip05+dFdpar_i)*hi;
            rows = rows+OCBVP.numode(arc);   % next equation

            Ji = Jip1;
            nrmJi = nrmJip1;
            dFdpar_i = dFdpar_ip1;
            nrmdFdpar_i = nrmdFdpar_ip1;
        end
    else % use analytical Jacobian

        [Ji,dFdpar_i] = odejac(xreg(1),yreg(:,1),FcnArgs{:});

        %%%%%%%
        %test
        %%%%%
%                 OCBVP.Fref=Freg(:,1);
%                 [Jitest,nFcalls,dFdpar_itest]=Fnumjac(ode,{xreg(1),yreg(:,1),FcnArgs{:}});
%                 if any(abs(Jitest(:)-Ji(:))/max(abs(Ji(:)))>1e0)
%                     %if any(abs(Jip05test(:)-Jip05(:))>1e-6)
%                     Ji
%                 end
%                 if any(abs(dFdpar_i(:)-dFdpar_itest(:))/max(abs(dFdpar_i(:)))>1e0)
%                     %if any(abs(dFdpar_ip05(:)-dFdpar_ip05test(:))>1e-6)
%                     dFdpar_i
%                  end
% 
        for i = 1:OCBVP.Nint(arc)
            % the left mesh point
            xi = xreg(i);
            yi = yreg(:,i);
            Fi = Freg(:,i);
            % the right mesh point
            xip1 = xreg(i+1);
            yip1 = yreg(:,i+1);
            Fip1 = Freg(:,i+1);
            [Jip1, dFdpar_ip1] = odejac(xip1,yip1,FcnArgs{:});
            % the midpoint
            hi = hreg(i);
            xip05 = (xi + xip1)/2;
            yip05 = (yi + yip1)/2 - hi/8*(Fip1 - Fi);
            [Jip05, dFdpar_ip05] = odejac(xip05,yip05,FcnArgs{:});  % recompute the Jacobian
            % test
            if OCMATCONT.OPTIONS.testmodus
                OCBVP.Fref=Fip1;
                [Jip1test,nFcalls,dFdpar_ip1test]=Fnumjac(ode,{xip1,yip1,FcnArgs{:}});
                if any(abs(Jip1test(:)-Jip1(:))/max(abs(Jip1(:)))>1e-6) %any(abs(Jip1test(:)-Jip1(:))>1e-5) %
                    %if any(abs(Jip05test(:)-Jip05(:))>1e-6)
                        Jip1
                end
                if any(abs(dFdpar_ip1(:)-dFdpar_ip1test(:))/max(abs(dFdpar_ip1(:)))>1e-7) %any(abs(dFdpar_ip1(:)-dFdpar_ip1test(:))>1e-5) %
                    %if any(abs(dFdpar_ip05(:)-dFdpar_ip05test(:))>1e-6)
                    dFdpar_ip1
                end
                OCBVP.Fref=Fmidreg(:,i);
                numJacOpt.diffvar=2;
                %             numJacOpt.vectvars=[];
                %             Jip05test=numjaccsd(ode,{xip05,yip05,FcnArgs{:}},length(yip05),numJacOpt);
                %             numJacOpt.diffvar=4;
                %             dFdpar_ip05test=numjaccsd(ode,{xip05,yip05,FcnArgs{:}},length(yip05),numJacOpt);
                [Jip05test,nFcalls,dFdpar_ip05test]=Fnumjac(ode,{xip05,yip05,FcnArgs{:}});
                if any(abs(Jip05test(:)-Jip05(:))/max(abs(Jip05(:)))>1e-1) %any(abs(Jip05test(:)-Jip05(:))>1e-5) %
                    %if any(abs(Jip05test(:)-Jip05(:))>1e-6)
                    Jip05
                end
                if any(abs(dFdpar_ip05(:)-dFdpar_ip05test(:))/max(abs(dFdpar_ip05(:)))>1e-7) %any(abs(dFdpar_ip05(:)-dFdpar_ip05test(:))>1e-5) %
                    %if any(abs(dFdpar_ip05(:)-dFdpar_ip05test(:))>1e-6)
                    dFdpar_ip05
                end

            end

            twiceJip05 = 2*Jip05;
            % assembly
            Jac(rows,cols) = -(OCBVP.In{arc}+hi/6*(Ji+twiceJip05*(OCBVP.In{arc}+hi/4*Ji)));
            cols = cols + OCBVP.numode(arc);
            Jac(rows,cols) = OCBVP.In{arc}-hi/6*(Jip1+twiceJip05*(OCBVP.In{arc}-hi/4*Jip1));
            %             last_cols(rows,:) = -hi*dFdpar_ip05 + hi^2/12*Jip05* ...
            %                 (dFdpar_ip1-dFdpar_i);
            last_cols(rows,:) = -1/6*(dFdpar_ip1-1/2*Jip05*(dFdpar_ip1-dFdpar_i)*hi+4*dFdpar_ip05+dFdpar_i)*hi;
            rows = rows+OCBVP.numode(arc);   % next equation

            Ji = Jip1;
            dFdpar_i = dFdpar_ip1;
        end
    end

    % Right BC
    Jac(1:OCBVP.nBCs,cols) = dGdyb(:,sum(OCBVP.numode(1:arc-1))+(1:OCBVP.numode(arc)));
    if arc<OCBVP.numarc
        cols = cols(end) + (1:OCBVP.numode(arc+1));
        rows = rows-OCBVP.numode(arc);
        rows = rows(end) + (1:OCBVP.numode(arc+1));
    end
end
if OCMATCONT.OPTIONS.messagedisplay
    ocmatmsg('elapsed time (J)=%.1f secs\n', etime(clock, StartTime));
end

Jac(:,end-OCBVP.npar+1:end) = last_cols;  % accumulated
