function [tcolmesh,coeff,tangent,iter]=newtcorr4dae(tcolmesh0,coeff0,tangent0)
% Newton solver and mesh adaptation, based on:
% Matlab Code BVPSUITE2.0 (beta-version July 2018).
% Authors: Winfried Auzinger, Merlin Fallahpour, Georg Kitzhofer,
% Othmar Koch, Gernot Pulverer, Gustaf Söderlind, Ewa B. Weinmueller and
% Stefan Wurm.
% Institute for Analysis and Scientific Computing, Vienna University of
% Technology, Vienna, Austria.


global OCMATCONT

% coefficientindex : last index of the coefficients excluding the
% continuation parameter
if isempty(tangent0)
    coefficientindex=1:OCMATCONT.MESHDATA.continuationindex-1;
else
    coefficientindex=1:OCMATCONT.MESHDATA.continuationindex;

    EN=[];
    EN(OCMATCONT.MESHDATA.continuationindex,1)=1;
end
coeff=coeff0;
tcolmesh=tcolmesh0; % time mesh cannot be changed during these steps, but can be set empty if no solution is detected

fcount = 0;   % number of function evaluations
jcount=1; % number of jacobian evaluations
iter= 1;   % number of iterations
lambda=1;
if OCMATCONT.OPTIONS.contlog
    jtime=0;
    ftime=0;
    showvariable={'method','iter','fcount','jcount','coeff_dist0','coeff_dist','tolfactor','lambda','l1','l2','ncorrections'};

    st=dbstack('-completenames');
    lst=length(st);
    OCMATCONT.writelogfile(sprintf('\nNewton solver: %s is called from function : ',st(1).name));
    if lst>=3
        for jj=lst:-1:3
            OCMATCONT.writelogfile(sprintf('%s -> ',st(jj).name));
        end
    end
    tmpshowv=showvariable;
    for ii=1:length(showvariable)
        if ii<length(showvariable)
            tmpshowv{ii}=[tmpshowv{ii} '\t'];
        else
            tmpshowv{ii}=[tmpshowv{ii} '\n'];
        end
    end
    OCMATCONT.writelogfile(sprintf(['\n' [tmpshowv{:}]]));
end

% depending on the convergence behavior of the initial data the initial
% Newton method is determined (method: 1:fast frozen or 2:line search, see
% schoebinger2015)
[coeff_new,tangent1,method,tolfactor,U,L,F,coeff_del,coeff_dist0,simplified_coeff_del,coeff_dist,fcount]= ...
    determineinitialmethod(tcolmesh0,coeff0,tangent0,fcount);

if OCMATCONT.OPTIONS.contlog
    OCMATCONT.writelogfile(sprintf('\nFunction: %s returns : \n','determineinitialmethod'));
    existidx=ismember(showvariable,who);
    logstr=writelogdata(showvariable,existidx);
    OCMATCONT.writelogfile(sprintf([logstr '\n']));
end

tangent=tangent1;

if tolfactor<1
    coeff(coefficientindex,1)=coeff_new(coefficientindex,1);
    return
end
if method==1
    coeff(coefficientindex,1)=coeff_new(coefficientindex,1);
end

while 1
    % THE MAIN LOOP:
    if iter>OCMATCONT.OPTIONS.maxnewtiter
        tcolmesh=[];
        coeff=[];
        tangent=[];
        return
    end

    if fcount>OCMATCONT.OPTIONS.maxprobes
        tcolmesh=[];
        coeff=[];
        tangent=[];
        return
    end
    switch method
        case 1
            needGlobJac=coeff_dist>OCMATCONT.OPTIONS.updatejacfactor*coeff_dist0;
            
            if needGlobJac
                DF=OCMATCONT.frechetder(tcolmesh,coeff,tangent,OCMATCONT.dae,OCMATCONT.bc,OCMATCONT.icfun,OCMATCONT.daejac,OCMATCONT.bcjac,OCMATCONT.icfunjac);

                [L,U]=lu(DF);            % LU-decomposition of DF
                coeff_del = U\(L\(- F));
                coeff_dist0=norm(coeff_del);
            else
                coeff_del=simplified_coeff_del;
                coeff_dist0=coeff_dist;
            end
            coeff_new(coefficientindex,1)=coeff(coefficientindex,1)+coeff_del(coefficientindex,1);
            new_F=OCMATCONT.operatoreq(tcolmesh,coeff_new,tangent,OCMATCONT.dae,OCMATCONT.bc,OCMATCONT.icfun);

            if OCMATCONT.OPTIONS.continuationmethod==1 && ~isempty(tangent)% Moore Penrose
                % find the Newton direction
                D=U\(L\[-new_F EN]);
                simplified_coeff_del=D(:,1);
                tangent=D(:,2);
                tangent=tangent/norm(tangent);
            else % Arclength
                simplified_coeff_del = U\(L\(- new_F));
            end

            coeff_dist=norm(simplified_coeff_del);

            if coeff_dist<(1-OCMATCONT.OPTIONS.lambdamin/2)*coeff_dist0 % Approximation improved
                tolfactor=check_tolerances(coeff_new,simplified_coeff_del);
                coeff(coefficientindex,1)=coeff_new(coefficientindex,1);         % Accept new approximation
                F=new_F;
                if OCMATCONT.OPTIONS.contlog
                    existidx=ismember(showvariable,who);
                    logstr=writelogdata(showvariable,existidx);
                    OCMATCONT.writelogfile(sprintf([logstr '\n']));
                end
                if tolfactor<1
                    return
                end
            elseif needGlobJac==0
                % try it again with new Jacobian
                method=1;
                if OCMATCONT.OPTIONS.contlog
                    existidx=ismember(showvariable,who);
                    logstr=writelogdata(showvariable,existidx);
                    OCMATCONT.writelogfile(sprintf([logstr '\n']));
                end
            else
                method=2;
                lambda=1;
                if OCMATCONT.OPTIONS.contlog
                    existidx=ismember(showvariable,who);
                    logstr=writelogdata(showvariable,existidx);
                    OCMATCONT.writelogfile(sprintf([logstr '\n']));
                end
            end
            iter=iter+1;

        case 2
            accept_lambda=coeff_dist<(1-OCMATCONT.OPTIONS.lambdamin/2)*coeff_dist0;
            ncorrections = 0;
            while ~accept_lambda
                if ncorrections == 0 % First correction, quadratical method
                    % We know G(0), G'(0) and G=G(lambda) by previous computations.
                    % To obtain a correction, we interpolate a quadratic polynomial through
                    % these values and define the minimum of this polynomial as the corrected lambda
                    % The respective polynomial reads G(x)=G(0) + x * G'(0) + x^2/lambda^2 *
                    % (G(lambda) - G(0) - lambda * G'(0)). Zeroing the first derivative yields:

                    % store lambda (=lambda_pred) for possible further correction
                    l2 = lambda;

                    coeff_distprime0 = -coeff_dist0;
                    lambda = - lambda^2 * coeff_distprime0 / (2* (coeff_dist - coeff_dist0 - lambda * coeff_distprime0));

                    % Allow lambda_cor between lambda_pred / 10 and 1
                    lambda = max(l2/10 , min(lambda, 1));

                    %***** Prepare for possible further correction (cubical, match notation with Num. Recipes)
                    l1 = lambda;

                    ncorrections = 1;
                    if OCMATCONT.OPTIONS.contlog
                        existidx=ismember(showvariable,who);
                        logstr=writelogdata(showvariable,existidx);
                        OCMATCONT.writelogfile(sprintf([logstr '\n']));
                    end
                else % Subsequent corrections, cubical method
                    cubcoeff = 1/(l1-l2) * [1/l1^2 -1/l2^2 ; -l2/l1^2 l1/l2^2] * ...
                        [coeff_dist1 - coeff_distprime0*l1 - coeff_dist0 ; coeff_dist2 - coeff_distprime0*l2 - coeff_dist0];

                    ca = cubcoeff(1);
                    cb = cubcoeff(2);

                    lambda = (-cb + sqrt(cb^2-3*ca*coeff_distprime0))/(3*ca);

                    % Allow lambda_cor between lambda_cor_old / 10 and lambda_cor_old /2
                    lambda = max(l1/10 , min(lambda, l1/2));

                    % Prepare for possible further correction
                    l2 = l1;
                    l1 = lambda;

                    ncorrections=ncorrections + 1;
                    if OCMATCONT.OPTIONS.contlog
                        existidx=ismember(showvariable,who);
                        logstr=writelogdata(showvariable,existidx);
                        OCMATCONT.writelogfile(sprintf([logstr '\n']));
                    end
                end

                if lambda<OCMATCONT.OPTIONS.lambdamin
                    method=3;  % switch to Trust Region Method
                    break;
                end

                % coeff_del does not change within this loop
                coeff_new(coefficientindex,1)=coeff(coefficientindex,1)+lambda*coeff_del(coefficientindex,1);  % get new approximation

                F=OCMATCONT.operatoreq(tcolmesh,coeff_new,tangent,OCMATCONT.dae,OCMATCONT.bc,OCMATCONT.icfun);
                fcount = fcount + 1;
                simplified_coeff_del = U\(L\-F);
                coeff_dist_old = coeff_dist; % store old coeff_dist for possible further correction
                coeff_dist=norm(simplified_coeff_del);

                accept_lambda=coeff_dist<(1-OCMATCONT.OPTIONS.lambdamin/2)*coeff_dist0;

                % ***** Prepare for possible further correction
                coeff_dist1 = coeff_dist;
                coeff_dist2 = coeff_dist_old;
                if OCMATCONT.OPTIONS.contlog
                    existidx=ismember(showvariable,who);
                    logstr=writelogdata(showvariable,existidx);
                    OCMATCONT.writelogfile(sprintf([logstr '\n']));
                end
            end

            if accept_lambda % a lambda has been accepted
                coeff(coefficientindex,1)=coeff_new(coefficientindex,1);    % Accept new solution

                if coeff_dist < OCMATCONT.OPTIONS.switchtoffnfactor * coeff_dist0
                    method= 1;
                else
                    method= 2;
                end
            end

            iter=iter+1;
            % At this point, method has been determined. Check Termination conditions if method~= 3
            if method~=3
                tolfactor=check_tolerances(coeff,simplified_coeff_del);
                if OCMATCONT.OPTIONS.contlog
                    existidx=ismember(showvariable,who);
                    logstr=writelogdata(showvariable,existidx);
                    OCMATCONT.writelogfile(sprintf([logstr '\n']));
                end
                if tolfactor < 1
                    coeff(coefficientindex,1)=coeff(coefficientindex,1)+simplified_coeff_del(coefficientindex,1);
                    if ~isempty(tangent)
                        DF=OCMATCONT.frechetder(tcolmesh,coeff,EN,OCMATCONT.dae,OCMATCONT.bc,OCMATCONT.icfun,OCMATCONT.daejac,OCMATCONT.bcjac,OCMATCONT.icfunjac);
                        %DF(OCMATCONT.MESHDATA.continuationindex,:)=tangent(:).';
                        [L,U]=lu(DF);            % LU-decomposition of DF


                        % find the Newton direction
                        tangent=U\(L\EN);
                        tangent=tangent/norm(tangent);
                    end
                    return
                end
            end


            % If we keep up doing the Line Search, we have to predict lambda for the next step
            if method== 2

                % ********** Keep old values of coeff_del for prediction
                coeff_dist0_old = coeff_dist0;

                % ********** Determine coeff_del
                DF=OCMATCONT.frechetder(tcolmesh,coeff,tangent,OCMATCONT.dae,OCMATCONT.bc,OCMATCONT.icfun,OCMATCONT.daejac,OCMATCONT.bcjac,OCMATCONT.icfunjac);
                [L,U]=lu(DF);                  % LU-decomposition of DF
                jcount=jcount+1;
                if OCMATCONT.OPTIONS.continuationmethod==1 && ~isempty(tangent)% Moore Penrose
                    % find the Newton direction
                    D = U\(L\[-F EN]);         % Newton correction
                    coeff_del=D(:,1);
                    tangent=D(:,2);
                    tangent=tangent/norm(tangent);
                else % Arclength
                    coeff_del = U\(L\(-F));
                end

                coeff_dist0 = norm(coeff_del);

                % ********** Predict lambda [Deuflhardt74]
                %dpr%      warning off MATLAB:divideByZero  % If Jacobians are the same, coeff_del = simplified_coeff_del
                mu = lambda * coeff_dist0_old / norm(simplified_coeff_del - coeff_del);
                %dpr%      warning on MATLAB:divideByZero

                if mu > 0.7
                    lambda = 1;
                else
                    lambda = mu;
                end

                % ********** Make sure all necessary variables are available for correction step
                coeff_new(coefficientindex,1)=coeff(coefficientindex,1)+lambda*coeff_del(coefficientindex,1);

                F=OCMATCONT.operatoreq(tcolmesh,coeff_new,tangent,OCMATCONT.dae,OCMATCONT.bc,OCMATCONT.icfun);
                fcount = fcount + 1;
                
                simplified_coeff_del = U\(L\-F);

                coeff_dist = norm(simplified_coeff_del);
                if OCMATCONT.OPTIONS.contlog
                    existidx=ismember(showvariable,who);
                    logstr=writelogdata(showvariable,existidx);
                    OCMATCONT.writelogfile(sprintf([logstr '\n']));
                end
            end
        case 3
            tcolmesh=[];
            coeff=[];
            tangent=[];
            if OCMATCONT.OPTIONS.contlog
                st=dbstack('-completenames');
                OCMATCONT.writelogfile(sprintf('\nNewton solver failed. Called from function %s at line %d\n',st(3).name,st(3).line));
            end
            return
    end
end

function [coeff_new,new_tangent,method,tolfactor,U,L,F,coeff_del0,coeff_dist0,coeff_del1,coeff_dist1,fcount]=determineinitialmethod(tcolmesh,coeff,tangent,fcount)
global OCMATCONT

if isempty(tangent) % if tangent is empty the continuaton parameter is fixed
    coefficientindex=1:OCMATCONT.MESHDATA.continuationindex-1;
    new_tangent=[];
else
    coefficientindex=1:OCMATCONT.MESHDATA.continuationindex;

    EN=[];
    EN(OCMATCONT.MESHDATA.continuationindex,1)=1;
end
coeff_new=coeff;

% first step: calculate Jacobian and residual
DF=OCMATCONT.frechetder(tcolmesh,coeff,tangent,OCMATCONT.dae,OCMATCONT.bc,OCMATCONT.icfun,OCMATCONT.daejac,OCMATCONT.bcjac,OCMATCONT.icfunjac);

F=OCMATCONT.operatoreq(tcolmesh,coeff,tangent,OCMATCONT.dae,OCMATCONT.bc,OCMATCONT.icfun);
fcount=fcount+1;

[L,U]=lu(DF);                  % LU-decomposition of DF
if ~isempty(tangent)
    D=U\(L\[-F EN]);
    coeff_del0=D(:,1);
    new_tangent=D(:,2);
    new_tangent=new_tangent/norm(new_tangent);
else
    coeff_del0=U\(L\(-F));
end

coeff_dist0=norm(coeff_del0);            % Norm of the 1. Newton correction

coeff_new(coefficientindex,1)=coeff(coefficientindex,1)+coeff_del0(coefficientindex,1); % full Newton step

tolfactor=check_tolerances(coeff_new,coeff_del0);
if tolfactor < 1
    method=1;  % Fast Frozen Newton
    coeff_del1=[];
    coeff_dist1=[];
    return
end
F=OCMATCONT.operatoreq(tcolmesh,coeff_new,tangent,OCMATCONT.dae,OCMATCONT.bc,OCMATCONT.icfun);   % new residual

fcount=fcount+1;

% calculate new deviation vector with new residual and fixed Jacobian
coeff_del1=U\(L\(-F));
coeff_dist1=norm(coeff_del1);

if coeff_dist1<(1-OCMATCONT.OPTIONS.lambdamin/2)*coeff_dist0 % Approximation improved
    if coeff_dist1<OCMATCONT.OPTIONS.switchtoffnfactor*coeff_dist0         % Approximation improved significantly
        method=1;             % -> Fast Frozen Newton

        if OCMATCONT.OPTIONS.contlog
            st=dbstack('-completenames');
            OCMATCONT.writelogfile(sprintf('%s : Solution improved significantly\n',st(1).name));
        end
        % See if tolerances are satisfied
        tolfactor=check_tolerances(coeff_new,coeff_del1);
        if tolfactor<1
            % Perform one last step towards the solution, since
            % approximatioin improved significantly take full Newtonstep
            % for new coefficient
            coeff_new(coefficientindex,1)=coeff_new(coefficientindex,1)+coeff_del1(coefficientindex,1);
            return
        end
    else % Approximation improved slightly
        method=2;             % -> Predictor Corrector Line Search
        if OCMATCONT.OPTIONS.contlog
            st=dbstack('-completenames');
            OCMATCONT.writelogfile(sprintf('%s : Solution improved slightly\n',st(1).name));
        end
    end
else % try another Newton step with smallest possible lambda

    coeff_new(coefficientindex,1)=coeff(coefficientindex,1)+OCMATCONT.OPTIONS.lambdamin*coeff_del0(coefficientindex,1);      % new approximation

    Fmin=OCMATCONT.operatoreq(tcolmesh,coeff_new,tangent,OCMATCONT.dae,OCMATCONT.bc,OCMATCONT.icfun);   % new residual
    fcount=fcount+1;

    coeff_del1_min = U\(L\(-Fmin));
    coeff_dist1min = norm(coeff_del1_min);

    if coeff_dist1min<(1-OCMATCONT.OPTIONS.lambdamin/2)*coeff_dist0 % Approximation improved
        method=2; % Predictor Corrector Line Search
    else
        method=3; % -> Trust Region Method
    end
    if OCMATCONT.OPTIONS.contlog
        st=dbstack('-completenames');
        OCMATCONT.writelogfile(sprintf('%s : Solution accepted for lambdamin\n',st(1).name));
    end
end

function tolfactor=check_tolerances(coeff,coeff_del)
global OCMATCONT
% Calculates the factor by which delta has to be reduced to satisfy the
% tolerances.
% out <   1   => Tolerances satisfied
% out >=  1   => Tolerances not satisfied

infnorm_x=max(abs(coeff));

infnorm_delta=max(abs(coeff_del));

% ***** check if tolerances are satisfied
tolfactor = max(infnorm_delta./(OCMATCONT.OPTIONS.newtonabstol+infnorm_x*OCMATCONT.OPTIONS.newtonreltol));

