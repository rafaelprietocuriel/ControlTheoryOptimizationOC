function [tmesh1,coeff,tangent,meshadapted]=meshadaptation_sbvpoc(tmesh1,coeff,tangent)


% meshadaptation ... calculate solution using adaptive meshes
%
% meshadaptation2(aTOL,rTOL,K,bvpfile,plotsol,tmesh1,coeff,bvpopt,plotres,numcols
% ax
% iter) calculates the solution of a BVP, given in
% bvpfile using an adaptive mesh technique in order to solve the prescribed
% tolerance requirements given in aTOL (absolute tolerance) and rTOL
% (relative tolerance).
%
% K limits the maximal ratio between the shortest and the longest
% subinterval. plotsol (0,1) toggles the plot of the solution after the
% solution process. An initial mesh and an initial solution can be passed
% to the routine using tmesh1 and coeff. If left empty, the mesh is determined
% from the bvpfile and the initial profile 0 is used. maxiter regulates the
% number of allowed adaptation steps.
%
% See also BVPSUITE, RUN
%
global OCMATCONT

aTOL=OCMATCONT.OPTIONS.meshadaptabstol;
rTOL=OCMATCONT.OPTIONS.meshadaptreltol;
K=OCMATCONT.OPTIONS.meshadaptk;
finemesh=OCMATCONT.OPTIONS.finemesh;
maxiter=OCMATCONT.OPTIONS.meshadaptmaxiter;


% Calculate mesh density rho from the mesh
rho=zeros(OCMATCONT.HE.TIMEDDATA.rightarcindex(end)-OCMATCONT.HE.numarc,1);
for arc=1:OCMATCONT.HE.numarc
    rho((OCMATCONT.HE.TIMEDDATA.leftarcindex(arc)-arc+1):(OCMATCONT.HE.TIMEDDATA.rightarcindex(arc)-arc))=1./diff(tmesh1(OCMATCONT.HE.TIMEDDATA.leftarcindex(arc):OCMATCONT.HE.TIMEDDATA.rightarcindex(arc)))/(OCMATCONT.HE.TIMEDDATA.rightarcindex(arc)-OCMATCONT.HE.TIMEDDATA.leftarcindex(arc))';
end
update_mode=1; % 1 .. update rho, 2 .. update N
meshadapted=0;

%M=OCMATCONT.HE.TIMEDDATA.nummeshintv;%length(tmesh1)-1; % number of intervals

M = OCMATCONT.HE.TIMEDDATA.nummeshintv;
N=M-1;
Ncompare = repmat(ceil(1e8/OCMATCONT.HE.numarc),OCMATCONT.HE.numarc,1); % initialized with a high value

%%%%
%%%%

% Variables
minimprove = 0.1; % 0.1 means at least 10 % points less to try further density adjustment
safetyN = 0; % 0.1 means at least 10 % more points than actually suggested - not accurate anymore - use safety_sigma instead
safety_sigma = 0.9; % safety factor used in the formula of N-preestimation
%LPFilter_k = order_method*2.75; % controls the gain of the adjustment of the mesh density - should be coupled to the method order
res_mode = 'tcolmid'; % triggers the choice of points for the evaluation of the residual

tmesh_intervalleft=tmesh1(OCMATCONT.HE.TIMEDDATA.leftarcindex);
tmesh_intervalright=tmesh1(OCMATCONT.HE.TIMEDDATA.rightarcindex);

maxnumeq=max([OCMATCONT.DOMAINDDATA(:).numeq]);
if length(aTOL) <= maxnumeq
    aTOL=repmat(aTOL(1),1,maxnumeq);
end
if length(rTOL) <= maxnumeq
    rTOL=repmat(rTOL(1),1,maxnumeq);
end

% DEFINITIONS ==================================
% Number of INTERVALS is M
% Number if internal GRID POINTS is N = M-1
% STEP SIZE on the auxiliary grid is dxi = 1/M
% rho is a vector with M components
% ==============================================

% Main loop: solution of problem + grid refinement

densityUpdates = 0;

for j=1:maxiter
    % ***************************
    % Calculate next mesh from rho
    %*****************************

    % *************
    % SOLVE PROBLEM
    % *************

    [coeff,tangent]=newtoncorrection(tmesh1,coeff,tangent);
    tau=tmesh1;
    
    [tcol ycol par]=discretizationdata(tmesh1,coeff,'collmeshdata');

    % indices for the solution at the actual grid
    leftarcindex1=OCMATCONT.HE.TIMEDDATA.leftarcindex;
    rightarcindex1=OCMATCONT.HE.TIMEDDATA.rightarcindex;
    nummeshintv1=OCMATCONT.HE.TIMEDDATA.nummeshintv;
    % **************
    % ESTIMATE ERROR
    % **************
    [ecol,maxecol,ecoeff,tcol2,ycol2,tau2,y2,tangent2,ecol2]=errorestimate(tmesh1,coeff,tangent,tcol,ycol,par);

    % indices for the solution at the finer grid
    
    % Calculate error norms - global error ;
    %ycol_norm = max(abs(ycol),[],1)
    %err_norm = max(abs(ecol),[],1)
    ycol_norm = abs(ycol);
    err_norm = abs(ecol);

    %ycol_norm2 = max(abs(ycol2),[],1);
    %err_norm2 = max(abs(ecol2),[],1);
    ycol_norm2 = abs(ycol2);
    err_norm2 = abs(ecol2);

    for tolt=1:size(err_norm,1);
        [qTol(tolt),jmax(tolt)] = max(err_norm(tolt,:) ./ (aTOL(tolt) + ycol_norm(tolt,:) .* rTOL(tolt))); % qTol is about tolfactor = errmax/tol
        [qTol2(tolt),jmax2(tolt)] = max(err_norm2(tolt,:) ./ (aTOL(tolt) + ycol_norm2(tolt,:) .* rTOL(tolt))); % qTol is about tolfactor = errmax/tol
    end

%     if plotres
%         disp(['Tolerance Factor: ' num2str(qTol)])
%     end
%     if finemesh
%         if plotres
%             disp(['Tolerance Factor on the fine grid: ' num2str(qTol2)])
%         end
%     end

    %%%%%%%%%%%%%%%%%%%%;
    % Update Procedure %;
    %%%%%%%%%%%%%%%%%%%%;
    
    %%%
    %%%
    %%
    % for testing purposes
    %%%
    %%%
    %%%
    %save
    
    if update_mode == 2; % update the number of mesh_points
        if qTol < 1
%             if plotres
%                 disp('Tolerance satisfied')
%             end
            tmesh1=tau;
            tolerated=0;
            break
        else
            if finemesh == 1
                if qTol2 < 1
%                     if plotres
%                         disp('Tolerance satisfied on the fine grid')
%                     end
                    coeff=ecoeff;
                    tmesh1=tau2;
                    tangent=tangent2;
                    tolerated=1;
                    meshadapted=2;
                    break
                end
            end
            if (update_mode==2 && densityUpdates <2)
                update_mode = 1;
                densityUpdates = 1;
            else
                for arc=1:OCMATCONT.HE.numarc
                    numcols=OCMATCONT.DOMAINDDATA(OCMATCONT.HE.arcindex(arc)).numcols;
                    N(arc) = ceil((1+safetyN)*N(arc)*(max(max(err_norm(leftarcindex1(arc):rightarcindex1(arc)),[],2)' ./ (safety_sigma*(aTOL))).^(1/(numcols+1))));
                end
            end
        end
    elseif update_mode == 1; % update the mesh density
        if qTol < 1
%             if plotres
%                 disp('Tolerance satisfied')
%             end
            tmesh1=tau;
            tolerated=0; % original time grid accepted
            break
        else
            if finemesh == 1
                if qTol2 < 1
%                     if plotres
%                         disp('Tolerance satisfied on the fine grid')
%                     end
                    coeff=ecoeff;
                    tmesh1=tau2;
                    tangent=tangent2;
                    meshadapted=2;
                    tolerated=1; % finer time grid accepted
                    break
                end
            end
            % Decide about the next step - Update N or rho ;
            Ncompare_old = Ncompare;
            counter=0;
            for arc=1:OCMATCONT.HE.numarc
                numcols=OCMATCONT.DOMAINDDATA(OCMATCONT.HE.arcindex(arc)).numcols;
                Ncompare(arc)=ceil((1+safetyN)*N(arc)*(max(max(err_norm(leftarcindex1(arc):rightarcindex1(arc)),[],2)' ./ (safety_sigma*(aTOL))).^(1/(numcols+1))));
                if update_mode == 1;
%                     if plotres
%                         disp(['N suggested: ' num2str(Ncompare)])
%                     end
                end
                if Ncompare(arc) > (1-minimprove)*Ncompare_old(arc) % check if density is finally adjusted
                    % here we have to do the following - it can happen, that
                    % the suggested number of points increases by a huge amount
                    % so it would be a bad idea to update the density in this
                    % last step - we have to use the older (better) density
                    % here - this may happen in rare cases
                    update_mode = 2;
                    if Ncompare(arc)<= Ncompare_old(arc)
                        N(arc) = Ncompare(arc);
                        Ncompare(arc) = 1e8;
                    else
                        N(arc) = Ncompare_old(arc);
                        Ncompare(arc)= 1e8;
                        rho((leftarcindex1(arc)-arc+1):(rightarcindex1(arc)-arc))=rhold((leftarcindex1(arc)-arc+1):(rightarcindex1(arc)-arc));
                        tmesh1(leftarcindex1(arc):rightarcindex1(arc))=x1old(leftarcindex1(arc):rightarcindex1(arc));
                    end
                else
                    if ~counter
                        densityUpdates = densityUpdates + 1;
                    end
                    counter=counter+1;
                end
            end
        end
    end

    % ******************
    % CALCULATE RESIDUAL
    % ******************
    if update_mode==1
        switch res_mode
            case 'mesh'
                resmesh = tmesh1;
                residual=abs(OCMATCONT.residual(tau,coeff,[],resmesh));
                % to get the residual for adjusting the mesh density we need a single
                % vector, so we take the 2-norm of the integral of it

                residualfinal=zeros(1,size(residual,2)-1);
                for i=1:size(residual,1)
                    residualadd=(residual(i,1:end-1)+residual(i,2:end)).*diff(tmesh1); % this means value on the lhs and rhs times the length
                    residualfinal=residualfinal+residualadd.^2;
                end
                residualfinal=sqrt(residualfinal);
            case 'tcol'
                resmesh = tcol;

                residual=abs(OCMATCONT.residual(tau,coeff,[],resmesh));
                % to get the residual for adjusting the mesh density we need a single
                % vector, so we take the 2-norm of the integral of it

                residualfinal=zeros(1,length(tmesh1)-1);
                for i=1:size(residual,1)
                    residualadd=(residual(i,1:end-1)+residual(i,2:end)).*diff(tcol);
                    for k=1:length(tmesh1)-1
                        residualadd2(k)=sum(residualadd(k*(order_method+1)-order_method:k*(order_method+1)));
                    end
                    residualfinal=residualfinal+residualadd2.^2;
                end
                residualfinal=sqrt(residualfinal);

            case 'midpoints'
                resmesh = tmesh1(1:end-1)+diff(tmesh1)/2;
                residual=abs(OCMATCONT.residual(tau,coeff,[],resmesh));
                residualfinal=zeros(1,size(residual,2));
                for i=1:size(residual,1)
                    residualadd=(residual(i,:)).*diff(tmesh1); % this is the value at the midpoint times the length
                    residualfinal=residualfinal+residualadd.^2;
                end
                residualfinal=sqrt(residualfinal);
            case 'tcolmid'
                resmesh = tcol(1:end-1)+diff(tcol)/2;
                [residual residualarcindex]=OCMATCONT.residual(tau,coeff,[],resmesh);
                residual=abs(residual);
                % to get the residual for adjusting the mesh density we need a single
                % vector, so we take the 2-norm of the integral of it
                arcindex=find(diff(tmesh1)==0);
                leftarcindex=[1 arcindex+1];
                rightarcindex=[arcindex numel(tmesh1)];
                arcindextcol=find(diff(tcol)==0);
                leftarcindextcol=[1 arcindextcol+1];
                rightarcindextcol=[arcindextcol numel(tcol)];
                difftcol=diff(tcol);
                residualfinal=zeros(1,length(tmesh1)-1);
                for arc=1:OCMATCONT.HE.numarc
                    arcindex=OCMATCONT.HE.arcindex(arc); % the arcindex of the actual arc
                    numcols=OCMATCONT.DOMAINDDATA(arcindex).numcols;
                    idx=residualarcindex(arc,1):residualarcindex(arc,2);
                    residualadd2=[];
                    for i=1:size(residual,1)
                        residualadd=residual(i,idx).*difftcol(idx);
                        counter=0;
                        for k=leftarcindex(arc):rightarcindex(arc)-1
                            counter=counter+1;
                            residualadd2(counter)=sum(residualadd(counter*(numcols+1)-numcols:counter*(numcols+1)))/(tmesh1(k+1)-tmesh1(k));
                        end
                        residualfinal((leftarcindex(arc)-arc+1):(rightarcindex(arc)-arc))=residualfinal(leftarcindex(arc):rightarcindex(arc)-1)+residualadd2.^2;
                    end
                    residualfinal=sqrt(residualfinal);
                end
        end
    end

    Nnew = N;
    Mnew = Nnew+1;

    % Display stuff
    rhold = rho;
    x1old = tmesh1;
    rhow_new=[];
    tmesh_new=[];
    coeff=[];
    for arc=1:OCMATCONT.HE.numarc
%         if plotres
%             disp(['N used:      ' num2str(Nnew(arc))])
%             disp(['---------------'])
%             disp(['Iteration ' num2str(j)])
%         end
        if update_mode == 1
%             if plotres
%                 disp(['Updating density'])
%             end
            fprintf(1,'%s%d\n','Density Update N=',N(arc));
            meshadapted=1;
        elseif update_mode == 2
%             if plotres
%                 disp('Updating N')
%             end
            fprintf(1,'%s%d\n','N Update N=',N(arc));
            meshadapted=2;
        end
        idx=(leftarcindex(arc)-arc+1):(rightarcindex(arc)-arc);
        idxtcol=(leftarcindextcol(arc)):(rightarcindextcol(arc));
        arcindex=OCMATCONT.HE.arcindex(arc); % the arcindex of the actual arc
        numcols=OCMATCONT.DOMAINDDATA(arcindex).numcols;
        LPFilter_k = numcols*2.75;
        if update_mode == 1;
            % Update density function
            rho(idx) = LPfilter(residualfinal(idx)',rho(idx),LPFilter_k);
        end;
        % check if the grid is within K-range
        maxrho = max(rho(idx));
        rho(idx)=max(rho(idx),maxrho/K);

        rho(idx) = bcTDFlogV4(rho(idx));           % Smoothing
        rho(idx) = rho(idx)*sum(1./rho(idx))/Mnew(arc);      % Normalize

        tmesh_arc=tmesh1(leftarcindex(arc):rightarcindex(arc));
        tmesh_arc = (tmesh_arc-tmesh_arc(1))/(tmesh_arc(end)-tmesh_arc(1));
        I = cumtrapz(tmesh_arc,([rho(idx(1)); rho(idx)]+[rho(idx); rho(idx(end))])'/2);
        tmesh_arc = pchip(I/I(end),tmesh_arc,0:1/(Mnew(arc)):1);
        rhow_new = [rhow_new;1./diff(tmesh_arc)'];
        tmesh_new= [tmesh_new tmesh_intervalleft(arc)+(tmesh_arc-0)*(tmesh_intervalright(arc)-tmesh_intervalleft(arc))/(1-0)]; %umrechnen von [0,1] auf [a,b]
        coeff=[coeff; ...
            initialmesh(tcol(idxtcol),ycol(:,idxtcol),tmesh_new,arc)];
    end
    coeff=[coeff;par(:)];
    rho=rhow_new;
    tmesh1=tmesh_new;
end % end j




if j==maxiter
    ocmaterror('The number of iterations allowed was succeeded - please check your settings and increase the maximum number of iterations');
    return
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

function [valerror,maxvalerror,coeff_2,tmesh1tau_2,val_2tmesh1undtau_2,tmesh1_2,val_2tmesh1_2,tangent_2,valerror2]=errorestimate(tmesh,coeff,tangent,collmesh,collval,par)
global OCMATCONT
domainddata=OCMATCONT.DOMAINDDATA;

leftarcindex1=OCMATCONT.HE.TIMEDDATA.leftarcindex;
rightarcindex1=OCMATCONT.HE.TIMEDDATA.rightarcindex;
OCMATCONT.meshadaptationflag=1;
nummesh=numel(tmesh);
% generate finer grid by adding midpoint in each grid interval
tmesh1_2(1:2:2*nummesh) = tmesh;
tmesh1_2(2:2:2*nummesh-1) = (tmesh(2:nummesh)+tmesh(1:nummesh-1))/2;
find(diff(tmesh1_2)==0);
tmesh1_2(ans(1:2:end))=[];
% for initial solution use interpolation of coarse solution on finer grid
diffmesh=diff(tmesh1_2);
arcindex=find(diffmesh==0);
leftarcindex=[1 arcindex+1];
rightarcindex=[arcindex numel(tmesh1_2)];
initcoeff=[];
counter=0;
for arc=1:OCMATCONT.HE.numarc
    %initcoeff=initialmesh(collmesh,collval,tmesh1_2,par);
    arcindex=OCMATCONT.HE.arcindex(arc);
    numcols=domainddata(arcindex).numcols;
    counter_start=counter+1;
    counter=counter+(numcols+1)*OCMATCONT.HE.TIMEDDATA.nummeshintv(arc)+1;
    % interpolate coefficients at hte finer time grid
    initcoeff=[initcoeff; ...
        initialmesh(collmesh(counter_start:counter),collval(:,counter_start:counter),tmesh1_2(leftarcindex(arc):rightarcindex(arc)),arc)];
end
initcoeff=[initcoeff;par(:)];
%calculate solution on the fine grid
[coeff_2,tangent_2]=newtoncorrection(tmesh1_2,initcoeff,tangent);
leftarcindex1_2=OCMATCONT.HE.TIMEDDATA.leftarcindex;
rightarcindex1_2=OCMATCONT.HE.TIMEDDATA.rightarcindex;

[tmesh1tau_2 val_2tmesh1undtau_2]=discretizationdata(tmesh1_2,coeff_2,'collmeshdata');
val_2tmesh1_2=discretizationdata(tmesh1_2,coeff_2,'meshdata');

% evaluate finer solution on coarser grid
%val_2tmesh1undtaut=calc_discretizationdata(tmesh1_2,coeff_2,'wert',1,collmesh); 

val_2tmesh1undtau=devalsbvpoc(coeff_2,tmesh1_2,collmesh); 
val_1tmesh1undtau_2=devalsbvpoc(coeff,tmesh,tmesh1tau_2); 

valerror=zeros(size(val_2tmesh1undtau));
valerror2=zeros(size(val_2tmesh1undtau_2));
for arc=1:OCMATCONT.HE.numarc
    arcindex=OCMATCONT.HE.arcindex(arc);
    numcols=domainddata(arcindex).numcols;
    valerror(:,leftarcindex1(arc):rightarcindex1(arc))=(2^numcols / (1-2^numcols)) * (val_2tmesh1undtau(:,leftarcindex1(arc):rightarcindex1(arc)) - collval(:,leftarcindex1(arc):rightarcindex1(arc)));
    % evaluate coarser solution on finer grid
    %val_1tmesh1undtau_2=calc_discretizationdata(tmesh,coeff,'value',1,tmesh1tau_2);

    valerror2(:,leftarcindex1_2(arc):rightarcindex1_2(arc))=(1 / (1-2^numcols)) * (val_2tmesh1undtau_2(:,leftarcindex1_2(arc):rightarcindex1_2(arc)) - val_1tmesh1undtau_2(:,leftarcindex1_2(arc):rightarcindex1_2(arc)));
end
% Definition der Maximumnorm der Absolutbeträge der Fehlervale
if length(valerror(:,1))>1
    maxvalerror=max(abs(valerror));
else
    maxvalerror=abs(valerror);
end

