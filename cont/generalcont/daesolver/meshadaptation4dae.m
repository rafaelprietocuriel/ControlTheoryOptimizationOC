function [tcolmeshnew,coeffnew,tangentnew,passed,maxerr_norm,maxerr_norm2]=meshadaptation4dae(tcolmesh,coeff,tangent)
global OCMATCONT
% meshadaptation ... calculate solution using adaptive meshes
%
% the meshadaptation algorithm is that from bvpsuite2
% meshadaptation2(aTOL,rTOL,K,bvpfile,plotsol,tmesh1,start,bvpopt,plotres,max
% iter) calculates the solution of a BVP, given in
% bvpfile using an adaptive mesh technique in order to solve the prescribed
% tolerance requirements given in aTOL (absolute tolerance) and rTOL
% (relative tolerance).
%
% K limits the maximal ratio between the shortest and the longest
% subinterval.
%

maxerr_norm=[];
maxerr_norm2=[];

sol=sol2daestruct(0,tcolmesh,coeff);

freeparameter=coeff(OCMATCONT.MESHDATA.freeparameterindex);
freeparameternum=OCMATCONT.freeparameternum;

update_mode=1; % 1 .. update rho, 2 .. update N

aTOL=OCMATCONT.OPTIONS.meshadaptabstol;
rTOL=OCMATCONT.OPTIONS.meshadaptreltol;
K=OCMATCONT.OPTIONS.krange;
maxiter=OCMATCONT.OPTIONS.meshadaptmaxiter;
finemesh=OCMATCONT.OPTIONS.meshadaptfinemesh;

% Variables
minimprove = 0.1; % 0.1 means at least 10 % points less to try further density adjustment
safetyN = 0; % 0.1 means at least 10 % more points than actually suggested - not accurate anymore - use safety_sigma instead
safety_sigma = 0.9; % safety factor used in the formula of N-preestimation
LPFilter_k = OCMATCONT.CollocationNumber*2.75; % controls the gain of the adjustment of the mesh density - should be coppled to the method order
res_mode='tcolmid'; % triggers the choice of points for the evaluation of the residual
M=length(sol.x)-1;
N=M-1;
N_orig=N;
Ncompare = 1e8; % initialised with a high value

tmesh1=sol.x;
left_end = tmesh1(1);
right_end = tmesh1(end);
h=diff(sol.x);
rho=(1./h/length(sol.x))';

if length(aTOL)<OCMATCONT.componentnumber
    aTOL=aTOL(ones(OCMATCONT.componentnumber,1));
end
if length(rTOL)<OCMATCONT.componentnumber
    rTOL=rTOL(ones(OCMATCONT.componentnumber,1));
end


% DEFINITIONS ==================================
% Number of INTERVALS is M
% Number of internal GRID POINTS is N = M-1
% STEP SIZE on the auxiliary grid is dxi = 1/M
% rho is a vector with M components
% ==============================================

% Main loop: solution of problem + grid refinement
if OCMATCONT.OPTIONS.contlog
    showvariable={'update_mode','iter','N_ori','N','Nnew','Ncompare','Ncompare_old','err_norm','err_norm2','qTol','qTol2','coeff_dist0','residualfinal'};
    tmpshowv=showvariable;
    for ii=1:length(showvariable)
        if ii<length(showvariable)
            tmpshowv{ii}=[tmpshowv{ii} '\t'];
        else
            tmpshowv{ii}=[tmpshowv{ii} '\n'];
        end
    end
end

densityUpdates = 0;
iter=0;
while 1
    iter=iter+1;
    if iter>maxiter
        break
    end
    % ***************************
    % Calculate next mesh from rho
    %*****************************
    if iter>1
        % tmesh1 is the new mesh,
        tcolmesh=makecollocationmesh(tmesh1);
        coeff=points2coeff(sol,tmesh1);
        % since the continuation parameter is fixed it has to be assured
        % that it is appended on the new coeff.
        coeff(end-freeparameternum+1:end)=freeparameter;
        [tcolmesh1,coeff1]=newtoncorrection(tcolmesh,coeff,[]);
        if isempty(tcolmesh1)
            passed=0;
            tcolmeshnew=[];
            coeffnew=[];
            tangentnew=[];
            return
        end
        sol=sol2daestruct(0,tcolmesh1,coeff1);
        freeparameter=sol.parameters;
    end
    % **************
    % ESTIMATE ERROR
    % **************

    [sol_coarse,sol_fine]=errorestimate4dae(sol);
    if isempty(sol_coarse)
        passed=0;
        tcolmeshnew=[];
        coeffnew=[];
        tangentnew=[];
        return
    end
    ecol=sol_coarse.data.errest;
    %ecoeff=sol_coarse.data.coeff;
    tcol=sol_coarse.data.xcol;
    ycol=sol_coarse.data.ycol;
    tcolmesh1_2=sol_fine.data.xcol;
    ycol2=sol_fine.data.ycol;
    coeff1_2=sol_fine.data.coeff;
    tmesh1_2=sol_fine.x;
    %y2=sol_fine.y;
    ecol2=sol_fine.data.errest;
    
    % Calculate error norms - global error ;  Richardson Extrapolation,
    % shampinewatts1976, ascheretal1995 5.5.2, for a detailed analysis and
    % derivation of local, global errors see deuflhardbornemann2008.
    ycol_norm = abs(ycol);
    err_norm = abs(ecol);

    ycol_norm2 = abs(ycol2);
    err_norm2 = abs(ecol2);

    qTol = zeros(size(err_norm,1),1);
    jmax = qTol;
    qTol2 = qTol;
    jmax2 = qTol;
    for tolt=1:size(err_norm,1)
        [qTol(tolt),jmax(tolt)] = max(err_norm(tolt,:) ./ (aTOL(tolt) + ycol_norm(tolt,:) .* rTOL(tolt))); % qTol is about tolfactor = errmax/tol
        [qTol2(tolt),jmax2(tolt)] = max(err_norm2(tolt,:) ./ (aTOL(tolt) + ycol_norm2(tolt,:) .* rTOL(tolt))); % qTol is about tolfactor = errmax/tol
    end
    if OCMATCONT.OPTIONS.contlog
        existidx=ismember(showvariable,who);
        logstr=writelogdata(showvariable,existidx);
        OCMATCONT.writelogfile(sprintf(['\n' [tmpshowv{:}]]));
        OCMATCONT.writelogfile(sprintf([logstr '\n']));
    end
    if qTol < 1
        coeffnew=coeff;
        tcolmeshnew=tcolmesh;
        passed=1;
        maxerr_norm=max(err_norm(:));
        maxerr_norm2=max(err_norm2(:));
        if iter>1
            tangentnew=calculatenewtangent(tcolmeshnew,coeffnew,sign(tangent(end)));
        else
            tangentnew=tangent;
        end
        break
    end
    if finemesh == 1
        if qTol2 < 1
            coeffnew=coeff1_2;
            tcolmeshnew=tcolmesh1_2;
            updatedaecoefficients(tmesh1_2,'MESHDATA');
            updatedaecoefficients(makefinemesh(tmesh1_2),'MESHDATA1_2');
            tangentnew=calculatenewtangent(tcolmeshnew,coeffnew,sign(tangent(end)));
            maxerr_norm=max(err_norm(:));
            maxerr_norm2=max(err_norm2(:));
            passed=1;
            break
        end
    end

    %%%%%%%%%%%%%%%%%%%%;
    % Update Procedure %; pulvereretal2011
    %%%%%%%%%%%%%%%%%%%%;

    if update_mode == 2; % update the number of mesh_points
        if densityUpdates <2
            update_mode = 1;
            densityUpdates = 1;
        else
            %Nold = N;
            N=ceil((1+safetyN)*N*(max(max(err_norm,[],2)./(safety_sigma*(aTOL))).^(1/(OCMATCONT.CollocationNumber+1))));
            %N = max(N,ceil(1.5^(1/OCMATCONT.CollocationNumber)*Nold)); % max increase
            % of number of points
        end
    elseif update_mode == 1; % update the mesh density
        % Decide about the next step - Update N or rho ;
        Ncompare_old = Ncompare;
        Ncompare=ceil((1+safetyN)*N*(max(max(err_norm,[],2)./ (safety_sigma*(aTOL))).^(1/(OCMATCONT.CollocationNumber+1))));

        if Ncompare > (1-minimprove)*Ncompare_old % check if density is finally adjusted
            update_mode = 2;
            % here we have to do the following - it can happen, that
            % the suggested number of points increases by a huge amount
            % so it would be a bad idea to update the density in this
            % last step - we have to use the older (better) density
            % here - this may happen in rare cases
            if Ncompare <= Ncompare_old || (~exist('rhold','var'))
                N = Ncompare;
                Ncompare = 1e8;
            else
                N = Ncompare_old;
                Ncompare = 1e8;
                rho=rhold;
                tmesh1=tmesh1old;
            end
        else
            densityUpdates = densityUpdates + 1;
        end
    end

    % Only a maximal increase of mesh points is allowed
    if N>=N_orig*OCMATCONT.OPTIONS.meshadaptfactormax
        %fprintf('\n  Suggested number of mesh points: %i\n Number of mesh points allowed in this step: %i\n', N+2, max(N_orig,N_orig*meshFactorMax))
        tcolmeshnew=[];
        coeffnew=[];
        tangentnew=[];
        passed=0;
        if OCMATCONT.OPTIONS.contlog
            existidx=ismember(showvariable,who);
            logstr=writelogdata(showvariable,existidx);
            OCMATCONT.writelogfile(sprintf([logstr '\n']));
        end
        return
    end

    % ******************
    % CALCULATE RESIDUAL
    % ******************
    if update_mode==1
        switch res_mode
            case 'mesh'
                resmesh = tmesh1;
                residual=abs(OCMATCONT.residual(tcolmesh,coeff,[],resmesh,OCMATCONT.dae));
                h=diff(tmesh1);
                % to get the residual for adjusting the mesh density we need a single
                % vector, so we take the 2-norm of the integral of it
%                 residualfinal2=zeros(1,size(residual,2)-1);
%                 for i=1:size(residual,1)
%                     residualadd=(residual(i,1:end-1)+residual(i,2:end)).*diff(x1); % this means value on the lhs and rhs times the length
%                     residualfinal2=residualfinal2+residualadd.^2;
%                 end
%                 residualfinal2=sqrt(residualfinal2);

                residualfinal=sqrt(sum(((2*residual(:,1:end-1)+diff(residual,[],2)).*h(ones(OCMATCONT.componentnumber,1),:)).^2));
            case 'tcol'
                resmesh = tcol;

                residual=abs(OCMATCONT.residual(tcolmesh,coeff,[],resmesh,OCMATCONT.dae));
                %residual=abs(computeResidual(problem,coeff,par,tau,rhoColl,resmesh,lambda_p));
                % to get the residual for adjusting the mesh density we need a single
                % vector, so we take the 2-norm of the integral of it

                % at each component res_k, we calculate the integral
                % R_ik:=\int_{\tau_i}^{\tau_{i+1}} res_k(s) ds
                % R_i:=sqrt(sum_k R_ik^2)

                residualfinal=zeros((OCMATCONT.MESHDATA.meshNumber-1)*OCMATCONT.componentnumber,(OCMATCONT.CollocationNumber+1));
                tmp=zeros(OCMATCONT.componentnumber,OCMATCONT.CollocationNumber+1);
                residualadd=(2*residual(:,1:end-1)+diff(residual,[],2)).*diff(tcol(ones(OCMATCONT.componentnumber,1),:),[],2);
                residualadd=reshape(residualadd(:),OCMATCONT.componentnumber*(OCMATCONT.CollocationNumber+1),[]);
                for ii=1:OCMATCONT.MESHDATA.meshNumber-1
                    tmp(:)=residualadd(:,ii);
                    residualfinal((ii-1)*OCMATCONT.componentnumber+1:ii*OCMATCONT.componentnumber,:)=tmp;
                end
                residualfinal=sqrt(sum(reshape(sum(residualfinal,2),OCMATCONT.componentnumber,[]).^2));

            case 'midpoints'
                h=diff(tmesh1);
                resmesh = tmesh1(1:end-1)+h/2;
                residual=abs(OCMATCONT.residual(tcolmesh,coeff,[],resmesh,OCMATCONT.dae));
                residualfinal=sqrt(sum((residual.*h(ones(OCMATCONT.componentnumber,1),:)).^2));

            case 'tcolmid'
                h=diff(tmesh1);
                resmesh = tcol(1:end-1)+diff(tcol)/2;

                residual=abs(OCMATCONT.residual(tcolmesh,coeff,[],resmesh,OCMATCONT.dae));
                %residual=abs(computeResidual(problem,coeff,par,tau,rhoColl,resmesh,lambda_p));
                % to get the residual for adjusting the mesh density we need a single
                % vector, so we take the 2-norm of the integral of it

                % at each component res_k, we calculate the wheighted
                % integral
                % R_ik:=1/(\tau_{i+1}-\tau_i)\int_{\tau_i}^{\tau_{i+1}} res_k(s) ds
                % R_i:=sqrt(sum_k R_ik^2)

                residualfinal=zeros((OCMATCONT.MESHDATA.meshNumber-1)*OCMATCONT.componentnumber,(OCMATCONT.CollocationNumber+1));
                tmp=zeros(OCMATCONT.componentnumber,OCMATCONT.CollocationNumber+1);
                H=h(ones(OCMATCONT.componentnumber,1),:);
                residualadd=residual.*diff(tcol(ones(OCMATCONT.componentnumber,1),:),[],2);
                residualadd=reshape(residualadd(:),OCMATCONT.componentnumber*(OCMATCONT.CollocationNumber+1),[]);
                for ii=1:OCMATCONT.MESHDATA.meshNumber-1
                    tmp(:)=residualadd(:,ii);
                    residualfinal((ii-1)*OCMATCONT.componentnumber+1:ii*OCMATCONT.componentnumber,:)=tmp;
                end
                residualfinal=sqrt(sum(reshape(sum(residualfinal,2)./H(:),OCMATCONT.componentnumber,[]).^2));
        end
    end
    
    Nnew = N;
    Mnew = Nnew+1;

    rhold = rho;
    tmesh1old = tmesh1;
    if update_mode == 1;
        % Update density function
        rho = LPfilter(residualfinal',rho,LPFilter_k);
    end

    % check if the grid is within K-range
    maxrho = max(rho);
    rho=max(rho,maxrho/K);
    
    rho = bcTDFlogV4(rho);           % Smoothing
    rho = rho*sum(1./rho)/Mnew;      % Normalize
    
    tmesh1 = (tmesh1-tmesh1(1))/(tmesh1(end)-tmesh1(1));
    I = cumtrapz(tmesh1,([rho(1); rho]+[rho; rho(end)])'/2);
    tmesh1 = pchip(I/I(end),tmesh1,0:1/(Mnew):1);
    rho = 1./diff(tmesh1)';
    tmesh1 = left_end+(tmesh1-0)*(right_end-left_end)/(1-0); %umrechnen von [0,1] auf [a,b]
    if update_mode == 2
        updatedaecoefficients(tmesh1,'MESHDATA');
        tmesh11_2=makefinemesh(tmesh1);
        updatedaecoefficients(tmesh11_2,'MESHDATA1_2');
    end
    if OCMATCONT.OPTIONS.contlog
        existidx=ismember(showvariable,who);
        logstr=writelogdata(showvariable,existidx);
        OCMATCONT.writelogfile(sprintf([logstr '\n']));
    end

end



% **********************************************************
% **********************************************************
% Functions for processing the residual ********************
% **********************************************************
% **********************************************************

function newrho = LPfilter(err,rho,LPFilter_k)
% Generate step density profile update
% For use with 2pBVP solver for 1st order systems
% Euler-Lagrange optimal grid is generated using
% using deadbeat control law. Local error is
% equidistributed over the interval.

% err is already on staggered grid in this version
k = LPFilter_k;       % k = 3 for DBC, k = 4 conv filter (root 1/3)
M = length(rho);

% Process input error
scalederr = M*abs(err).^(1/k);
% scalederr = processV4(bcTDFlogV4(scalederr)); % SWITCHED OFF

ratio = scalederr;
% ratio is suggested update of rho; compute next rho
rhonew = rho.*ratio;
%rhonew = bcTDFlogV4(rhonew);                   % SWITCHED OFF
newrho = rhonew;%*sum(1./rhonew)/M;


function out = bcTDFlogV4(in)
% Boundary corrected Topelitz digital filter, multiplicative (logarithmic)
% Applies 2nd order convolution LP filter repeatedly
% to input signal until output signal is sufficiently smooth

smooth = 50;     % Define smoothness requirement
N = length(in);
old = in;
signal = old;    % Mem alloc
snratio = 0;

while snratio < smooth
    % Apply filter until S/N ratio at this stage is at least = smooth
    for i=2:N-1
        signal(i) = (old(i-1)*old(i)^2*old(i+1))^(1/4);
    end
    % Boundary corrections
    signal(1) = (old(1)^3*old(2)^3/old(3)^2/old(4)*old(5))^(1/4);
    % signal(1) = (old(1)^2*old(2)^3/old(3))^(1/4);
    signal(N) = (old(N)^3*old(N-1)^3/old(N-2)^2/old(N-3)*old(N-4))^(1/4);
    % signal(N) = (old(N)^2*old(N-1)^3/old(N-2))^(1/4);
    
    % Compute noise
    noise = old - signal;
    s = norm(signal);
    n = norm(noise);
    snratio = s/(n + 1e-3*s);
    old = signal;
end

out = signal;