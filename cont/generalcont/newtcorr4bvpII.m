function [coeff,tangent,itcount,logstruct,fcount]=newtcorr4bvpII(tmesh0,coeff0,tangent0)
%    newtcorr4bvpII(Fhandle,DFhandle,FDFhandle,x0,bvpfile,tau,bvpopt,parameters,psival,psi)
%   SOLVE_NONLINEAR_SYS  Solve nonlinear systems obtained by SBVPCOL and SBVPERR
%
%   This routine is private and should not be accessed by other callers than SBVPCOL and SBVPERR

% ********** store some fields of bvpopt as variables for quick reference
global OCMATCONT


AbsTol=OCMATCONT.OPTIONS.newtonabstol;
RelTol=OCMATCONT.OPTIONS.newtonreltol;
max_F_evals=10;%OCMATCONT.OPTIONS.maxprobes;
max_iter=OCMATCONT.OPTIONS.maxnewtiter;
display=OCMATCONT.OPTIONS.display;
log=OCMATCONT.OPTIONS.log;
TRM=OCMATCONT.OPTIONS.trm; % trust region method

% ********** some initializations
fcount = 1;   % number of function evaluations
itcount= 1;   % number of iterations
numdvariables=numel(coeff0);
new_coeff=zeros(numdvariables,1);
coeff=coeff0;
lambda=1;
lambdamin=OCMATCONT.OPTIONS.lambdamin;
nPreviousTRMIterates = 0;
logstruct = [];
R=[];
R(numdvariables,1)=1;

if log
    logstruct.DOC = [];
    logstruct.tol_factor = [];
    logstruct.G0 = [];
    logstruct.G = [];
    logstruct.jac_update = [];
    logstruct.lambda = [];
    logstruct.nCorrections = [];

    logstruct.runtime.DFeval = [];
    logstruct.runtime.Feval = [];
    logstruct.runtime.lu = [];
    logstruct.runtime.resubstitute = [];
end

% Structure of the Algorithm:
% The Zerofinder consists of 3 different algorithms:
%    *) Fast Frozen Newton (cheap, but small domain of convergence (DOC) G1)
%    *) Predictor-Corrector Linesearch (expensive, but large DOC G2)
%    *) Trust Region Method (even more expensive, but even larger DOC G3)
%
% For the choice of the appropriate algorithm to use, it is necessary to determine the
% position of the current approximation w.r.t. the DOCs G1, G2, G3

% ********** Determine position of initial approximation

if display=='i'  % display information every iteration step
    fprintf('\n Determining zerofinder algorithm ... \n\n');
end


[TolFactor,DOC,new_coeff,tangent,U,L,G0,G,F,delta_coeff,simplified_delta_coeff,simplified_tangent,fcount,logstruct] = ...
    determine_position(tmesh0,coeff0,tangent0,fcount,logstruct,R,[],[]);
%    determine_position(Fhandle,DFhandle,coeff,bvpfile,tau,bvpopt,parameters,psival,psi,lambdamin,fcount,logstruct,[],[]);

if TolFactor < 1
    if display ~= 'o'  % off
        fprintf('\n\n Tolerances satisfied\n');
    end
    return
end

if DOC == 1
    coeff = new_coeff; % Accept new approximation
end

% ********** initialize Zerofinder Display
if display=='i'  % display information every iteration step
    switch DOC
        case 1
            fprintf(' Fast Frozen Newton\n');
            fprintf(' STEP   F-COUNT   JAC_UPDATE   |DELTA_X|/|X|   IMP_FAC   TOLDIST_FAC \n');
        case 2
            fprintf(' Predictor-Corrector Line Search\n');
            fprintf(' STEP   F-COUNT     LAMBDA     |DELTA_X|/|X|   IMP_FAC   TOLDIST_FAC \n');
        case 3
            fprintf(' Trust Region Method\n');
    end
end

% ********** Choice-of-Algorithm Loop
while 1

    if log 
        logstruct.DOC = [logstruct.DOC DOC]; 
    end
    if itcount > max_iter
        error(' Maximum number of iterations exceeded');
    end
    if fcount > max_F_evals
        error(' Maximum number of function evaluations exceeded');
    end
    switch DOC

        % **************************************************************************
        case 1 % ******************* Fast Frozen Newton ****************************
            % **************************************************************************
            % variables that have to be available: a (=coeff), G0, G, F (=F(coeff)), simplified_delta_coeff

            % ********** Determine if Jacobian has to be updated
            UpdateJac=G>OCMATCONT.OPTIONS.updatejacfactor*G0;

            if log
                logstruct.lambda = [logstruct.lambda 1];
                logstruct.nCorrections = [logstruct.nCorrections 0];
                logstruct.jac_update = [logstruct.jac_update UpdateJac];
            end

            % ********** Determine Jacobian
            if UpdateJac % Update Jacobian
                cpt=cputime;
                DF=OCMATCONT.frechetder(tmesh0,coeff,tangent);
                DF=[DF;tangent(:).'];
                if log
                    logstruct.runtime.DFeval = [logstruct.runtime.DFeval cputime-cpt];
                    logstruct.condestDF = condest(DF);
                end
                cpt=cputime;
                [L,U]=lu(DF);                  % LU-decomposition of DF
                if log
                    logstruct.runtime.lu = [logstruct.runtime.lu cputime-cpt];
                end
                lastwarn('');                  % initialize warning state
                cpt=cputime;
                if OCMATCONT.OPTIONS.continuationmethod==1 % Moore Penrose
                    % find the Newton direction
                    D=U\(L\([-F R] ));
                    delta_coeff=D(:,1);
                    tangent=D(:,2);
                    tangent=tangent/norm(tangent);
                else % Arclength
                    delta_coeff = U\(L\(- F));         % Newton correction
                end
                if log logstruct.runtime.resubstitute = ...
                        [logstruct.runtime.resubstitute cputime-cpt];
                end

                if length(lastwarn)            % Exit if DF is singular
                    %      error(' System matrix is close to singular. Try a refined mesh or another initial approximation.');
                end

                G0 = norm(delta_coeff);            % Norm of the Newton correction
            else % Iterate with frozen Jacobian
                delta_coeff = simplified_delta_coeff;
                G0 = G; % = norm(simplified_delta_coeff)
            end

            % ********** Perform Newton step
            new_coeff=coeff+delta_coeff;      % new approximation

            cpt = cputime;
            new_F=OCMATCONT.operatoreq(tmesh0,new_coeff,tangent);   % new residual
            if log 
                logstruct.runtime.Feval = [logstruct.runtime.Feval cputime-cpt]; 
            end
            fcount = fcount +1;             % increase number of function evaluations
            % ********** Determine Position of new solution approximation
            cpt = cputime;
            if OCMATCONT.OPTIONS.continuationmethod==1 % Moore Penrose
                % find the Newton direction
                D=U\(L\([-F R] ));
                simplified_delta_coeff=D(:,1);
                simplified_tangent=D(:,2);
                simplified_tangent=simplified_tangent/norm(simplified_tangent);
            else % Arclength
                simplified_delta_coeff=U\(L\(- F));         % Newton correction
            end
            if log
                logstruct.runtime.resubstitute = [logstruct.runtime.resubstitute cputime-cpt];
            end
            G = norm(simplified_delta_coeff);
            if log
                logstruct.G0 = [logstruct.G0 G0];
                logstruct.G =  [logstruct.G G];
            end
            if G < (1-lambdamin/2) * G0 % Approximation improved
                % ***** Check if Tolerances are satisfied
                % ***** TolFactor: Factor by which the correction had to be
                % ***** reduced in order to satisfy the tolerances
                TolFactor = check_tolerances(new_coeff,simplified_delta_coeff);

                if log 
                    logstruct.tol_factor = [logstruct.tol_factor TolFactor]; 
                end

                if TolFactor < 1
                    % Perform the last step towards the solution
                    new_coeff=new_coeff+simplified_delta_coeff;

                    if display=='i'  % display information every iteration step
                        fprintf(' %3i   %5i    %8i     %13.2e    %.2e   %5.2e\n',...
                            itcount,fcount,UpdateJac,G/norm(coeff(:)),G0/G,TolFactor);
                    end

                    if display ~= 'o'  % off
                        fprintf('\n\n Tolerances satisfied\n');
                    end
                    DF=OCMATCONT.frechetder(tmesh0,coeff,tangent);
                    DF=[DF;simplified_tangent(:).'];
                    [L,U] = lu(DF);
                    tangent=U\(L\(R));
                    tangent=tangent/norm(tangent);
                    return
                end
                DOC = 1;           % Keep iterating with the Fast Frozen Newton
                coeff = new_coeff;         % Accept new approximation
                F = new_F;

                if display=='i'  % display information every iteration step
                    fprintf(' %3i   %5i    %8i     %13.2e    %.2e   %5.2e\n',...
                        itcount,fcount,UpdateJac,G/norm(coeff(:)),G0/G,TolFactor);
                end
            elseif UpdateJac == 0 % Jacobian has not been updated in the previous step
                DOC = 1;           % -> Try FFN once again,
                %    but update Jacobian this time
                %    (We would have to evaluate the Jacobian
                %    anyway in the Line Search Algorithm)

                if display=='i'  % display information every iteration step
                    fprintf(' %3i   %5i    %8i     %13.2e    %.2e   %5.2e   NOT ACCEPTED\n',...
                        itcount,fcount,UpdateJac,G/norm(coeff(:)),G0/G,TolFactor);
                end
            else
                DOC = 2;           % -> Switch to Predictor Corrector Line Search

                % ********** Prepare for next step with the respective method
                lambda = 1;

                % Use Jacobian computed here
                if log 
                    logstruct.jac_update  = [logstruct.jac_update  0]; 
                end

                if display=='i'
                    fprintf(' %3i   %5i    %8i     %13.2e    %.2e   %5.2e   NOT ACCEPTED\n',...
                        itcount,fcount,UpdateJac,G/norm(coeff(:)),G0/G,TolFactor);

                    fprintf('\n Switching to Predictor-Corrector Line Search\n');
                    fprintf(' STEP   F-COUNT     LAMBDA     |DELTA_X|/|X|   IMP_FAC   TOLDIST_FAC \n');
                end
            end

            itcount = itcount + 1;


            % **************************************************************************
        case 2 % *************** Predictor-Corrector Line Search *******************
            % **************************************************************************
            % Necessary variables at this point: a (=coeff), delta_coeff, simplified_delta_coeff,
            % lambda (=lambda_pred), G0, G (=G(lambda_pred))

            % Let G(lambda) = ||DF^(-1) * F(coeff + lambda*delta_coeff)||
            % At this point, coeff, delta_coeff, lambda, G(0), G'(0), G(lambda) are known
            % We check if lambda is accepted and correct lambda otherwise

            Accept_Lambda = G < (1-lambdamin/2)*G0;
            nCorrections = 0;

            while ~Accept_Lambda
                % ********** Repeatedly correct lambda until it is accepted or a termination condition is met

                if nCorrections == 0 % First correction, quadratical method
                    % We know G(0), G'(0) and G=G(lambda) by previous computations.
                    % To obtain a correction, we interpolate a quadratic polynomial through
                    % these values and define the minimum of this polynomial as the corrected lambda
                    % The respective polynomial reads G(coeff)=G(0) + coeff * G'(0) + coeff^2/lambda^2 *
                    % (G(lambda) - G(0) - lambda * G'(0)). Zeroing the first derivative yields:

                    % store lambda (=lambda_pred) for possible further correction
                    l2 = lambda;

                    Gprime0 = -G0;
                    lambda = - lambda * Gprime0 / (2* (G - G0 - lambda * Gprime0));

                    % Allow lambda_cor between lambda_pred / 10 and 1
                    lambda = max(l2/10 , min(lambda, 1));

                    %***** Prepare for possible further correction (cubical, match notation with Num. Recipes)
                    l1 = lambda;

                    nCorrections = 1;
                else % Subsequent corrections, cubical method
                    coeff = 1/(l1-l2) * [1/l1^2 -1/l2^2 ; -l2/l1^2 l1/l2^2] * ...
                        [G1 - Gprime0*l1 - G0 ; G2 - Gprime0*l2 - G0];

                    ca = coeff(1);
                    cb = coeff(2);

                    lambda = (-cb + sqrt(cb^2-3*ca*Gprime0))/(3*ca);

                    % Allow lambda_cor between lambda_cor_old / 10 and lambda_cor_old /2
                    lambda = max(l1/10 , min(lambda, l1/2));

                    % Prepare for possible further correction
                    l2 = l1;
                    l1 = lambda;

                    nCorrections = nCorrections + 1;
                end

                if lambda < lambdamin
                    DOC = 3;  % switch to Trust Region Method
                    break;
                end

                new_coeff=coeff+lambda*delta_coeff;  % get new approximation

                cpt = cputime;
                F=OCMATCONT.operatoreq(tmesh0,new_coeff,tangent);   % evaluate F at lambda_cor
                if log 
                    logstruct.runtime.Feval = [logstruct.runtime.Feval cputime-cpt]; 
                end

                fcount = fcount + 1;

                cpt = cputime;
                if OCMATCONT.OPTIONS.continuationmethod==1 % Moore Penrose
                    % find the Newton direction
                    D=U\(L\([-F R]));         % Newton correction
                    simplified_delta_coeff=D(:,1);
                    simplified_tangent=D(:,2);
                    simplified_tangent=simplified_tangent/norm(simplified_tangent);
                else % Arclength
                    simplified_delta_coeff = U\(L\(-F));            % simplified delta_coeff
                end
                if log 
                    logstruct.runtime.resubstitute=[logstruct.runtime.resubstitute cputime-cpt]; 
                end
                G_old = G; % store old G for possible further correction
                G = norm(simplified_delta_coeff);
                Accept_Lambda = G < (1-lambdamin/2) * G0;
                % ***** Prepare for possible further correction
                G1 = G;
                G2 = G_old;
            end

            if Accept_Lambda % a lambda has been accepted
                coeff = new_coeff;    % Accept new solution

                if G < OCMATCONT.OPTIONS.switchtoffnfactor * G0
                    DOC = 1;
                else
                    DOC = 2;
                end
            end

            itcount = itcount + 1;

            if log
                logstruct.lambda = [logstruct.lambda lambda];
                logstruct.nCorrections = [logstruct.nCorrections nCorrections];
                logstruct.G0 = [logstruct.G0 G0];
                logstruct.G =  [logstruct.G G];
            end

            % At this point, DOC has been determined. Check Termination conditions if DOC ~= 3
            if DOC ~= 3
                TolFactor = check_tolerances(coeff,simplified_delta_coeff);

                if log 
                    logstruct.tol_factor = [logstruct.tol_factor TolFactor]; 
                end

                if TolFactor < 1
                    new_coeff=coeff+simplified_delta_coeff;
                    tangent=simplified_tangent;
                    if display=='i'  % display information every iteration step
                        fprintf(' %3i   %5i        %5.4f    %12.2e    %.2e   %5.2e\n',...
                            itcount,fcount,lambda,G/norm(coeff(:)),G0/G,TolFactor);
                    end

                    if display ~= 'o'  % off
                        fprintf('\n\n Tolerances satisfied\n');
                    end

                    return
                end
            elseif log
                logstruct.tol_factor = [logstruct.tol_factor -1];
            end

            if display=='i'  % display information every iteration step
                fprintf(' %3i   %5i        %5.4f    %12.2e    %.2e   %5.2e\n',...
                    itcount,fcount,lambda,G/norm(coeff(:)),G0/G,TolFactor);
                if DOC == 1
                    fprintf('\n Switching to Fast Frozen Newton\n');
                    fprintf(' STEP   F-COUNT   JAC_UPDATE   |DELTA_X|/|X|   IMP_FAC   TOLDIST_FAC \n');
                elseif DOC == 3
                    fprintf('\n Switching to Trust Region Method\n');
                end
            end


            % If we keep up doing the Line Search, we have to predict lambda for the next step
            if DOC == 2

                % ********** Keep old values of delta_coeff for prediction
                G0_old = G0;

                % ********** Determine delta_coeff
                cpt = cputime;
                DF=OCMATCONT.frechetder(tmesh0,coeff,tangent);
                DF=[DF;tangent(:).'];
                if log
                    logstruct.runtime.DFeval = [logstruct.runtime.DFeval cputime-cpt];
                    logstruct.jac_update = [logstruct.jac_update 1];
                    logstruct.condestDF = condest(DF);
                end

                cpt=cputime;
                [L,U]=lu(DF);                  % LU-decomposition of DF
                if log 
                    logstruct.runtime.lu = [logstruct.runtime.lu cputime-cpt]; 
                end
                lastwarn('');                  % initialize warning state
                cpt=cputime;
                if OCMATCONT.OPTIONS.continuationmethod==1 % Moore Penrose
                    % find the Newton direction
                    D=U\(L\([-F R] ));
                    delta_coeff=D(:,1);
                    tangent=D(:,2);
                    tangent=tangent/norm(tangent);
                else % Arclength
                    delta_coeff = U\(L\(- F));         % Newton correction
                end
                if log logstruct.runtime.resubstitute = ...
                        [logstruct.runtime.resubstitute cputime-cpt]; end

                if length(lastwarn)            % Exit if DF is singular
                    %         error(' System matrix is close to singular. Try a refined mesh or another initial approximation.');
                end

                G0 = norm(delta_coeff);

                % ********** Predict lambda [Deuflhardt74]
                %dpr%      warning off MATLAB:divideByZero  % If Jacobians are the same, delta_coeff = simplified_delta_coeff
                mu = lambda * G0_old / norm(simplified_delta_coeff - delta_coeff);
                %dpr%      warning on MATLAB:divideByZero

                if mu > 0.7
                    lambda = 1;
                else
                    lambda = mu;
                end

                % ********** Make sure all necessary variables are available for correction step
                new_coeff=coeff+lambda*delta_coeff;

                cpt = cputime;
                F=OCMATCONT.operatoreq(tmesh0,new_coeff,tangent); 
                if log 
                    logstruct.runtime.Feval = [logstruct.runtime.Feval cputime-cpt]; 
                end

                fcount = fcount + 1;

                cpt = cputime;
                if OCMATCONT.OPTIONS.continuationmethod==1 % Moore Penrose
                    % find the Newton direction
                    D=U\(L\([-F R]));         % Newton correction
                    simplified_delta_coeff=D(:,1);
                    simplified_tangent=D(:,2);
                    simplified_tangent=simplified_tangent/norm(simplified_tangent);
                else % Arclength
                    simplified_delta_coeff = U\(L\(-F));            % simplified delta_coeff
                end
                if log 
                    logstruct.runtime.resubstitute=[logstruct.runtime.resubstitute cputime-cpt]; 
                end

                G = norm(simplified_delta_coeff);
            end


            % **************************************************************************
        case 3 % ********************* Trust Region Method *************************
            % **************************************************************************

            % Necessary variables at this point: a (=coeff), nPreviousTRMIterates
            %---------Achtung: Temporäre Änderung von mir!!!!
            if (~TRM)
                err('trm');
            end
            %Ende Änderung
            % ********** Determine accuracy for Trust Region Method
            switch nPreviousTRMIterates
                case 0
                    bvpopt.ZfOpt.TolX = sqrt(max(RelTol,AbsTol));
                    bvpopt.ZfOpt.TolFun = 0;        % Assemble necessary options for FSOLVE
                    bvpopt.ZfOpt.Jacobian = 'on';   % (Display, MaxIter, MaxFunEvals are user parameters)
                    bvpopt.ZfOpt.LargeScale = 'on';
                case 1
                    bvpopt.ZfOpt.TolX = max(RelTol,AbsTol);
                case 2
                    bvpopt.ZfOPt.TolX = 1000*eps;
                case 3
                    error(' Requested tolerance is beyond the accuracy of the algorithm');
            end

            bvpopt.ZfOpt.Display = 'iter';
            % *********** Limit the maximum number of iterations and FunEvals for TRM
            bvpopt.ZfOpt.MaxIter = max_iter - itcount;
            bvpopt.ZfOpt.MaxFunEvals = max_F_evals - fcount;

            % ********** Invoke Trust Region Method
            [coeff,dummy,F,ExitFlag,lsqLog,dummy2,DF] = ...
                lsqnonlin(FDFhandle,coeff,[],[],bvpopt.ZfOpt,bvpfile,tau,bvpopt,parameters,psival,psi);

            nPreviousTRMIterates = nPreviousTRMIterates +1;
            fcount = fcount + lsqLog.funcCount;
            icount = itcount + lsqLog.iterations;

            % ********** Intercept numerically singular Jacobians
            if length(findstr(lastwarn,'singular'))
                %      error(' System matrix is close to singular. Try a refined mesh or another initial approximation.');
            end

            % ***** Determine whether anything has gone wrong
            if ExitFlag <= 0
                error(' Zerofinder did not converge. Increasing MaxIter and MaxFunEvals might help.');
            end


            % ********** Determine position of new approximation
            [TolFactor, DOC, new_coeff, U, L, G0, G, F, delta_coeff, simplified_delta_coeff, fcount, logstruct] = ...
                determine_position(Fhandle,DFhandle,coeff,bvpfile,tau,bvpopt,parameters,psival,psi,lambdamin,fcount,logstruct,F,DF);

            if TolFactor < 1
                if display ~= 'o'  % off
                    fprintf('\n\n Tolerances satisfied\n');
                end
                return
            end

            % ***** We only log those quantities if we continue iterating
            if log
                logstruct.lambda = [logstruct.lambda -1];
                logstruct.nCorrections = [logstruct.nCorrections -1];
                logstruct.jac_update = [logstruct.jac_update -1];
                logstruct.G0 = [logstruct.G0 G0];
                logstruct.G =  [logstruct.G G];
                logstruct.tol_factor = [logstruct.tol_factor TolFactor];
            end


            if DOC == 1
                coeff = new_coeff; % Accept new approximation
            elseif DOC == 2  % use Jacobian computed here
                if log logstruct.jac_update  = [logstruct.jac_update  0]; end %dpr%
            end


            if display=='i'  % display information every iteration step
                switch DOC
                    case 1
                        fprintf(' Switching to Fast Frozen Newton\n\n');
                        fprintf(' STEP   F-COUNT   JAC_UPDATE   |DELTA_X|/|X|   IMP_FAC   TOLDIST_FAC \n');
                    case 2
                        fprintf(' Switching to Predictor-Corrector Line Search\n\n');
                        fprintf(' STEP   F-COUNT     LAMBDA     |DELTA_X|/|X|   IMP_FAC   TOLDIST_FAC \n');
                    case 3
                        fprintf(' Continuing with Trust Region Method\n\n');
                end
            end
    end
end



% *************************************************************************
% ********** DETERMINE POSITION OF INITIAL APPROXIMATION ******************
% *************************************************************************

function [TolFactor,DOC,new_coeff,tangent,U,L,G0,G,F,delta_coeff,simplified_delta_coeff,simplified_tangent,fcount,logstruct]= ...
    determine_position(tmesh,coeff,tangent,fcount,logstruct,R,F,DF)

global OCMATCONT
log=OCMATCONT.OPTIONS.log;

% ********** Determine Newton correction
% ***** Evaluate Jabobian, if not provided (by fsolve)
if isempty(F)
    cpt = cputime;
    DF=OCMATCONT.frechetder(tmesh,coeff,tangent);
    DF=[DF;tangent(:).'];
    if log
        logstruct.runtime.DFeval = [logstruct.runtime.DFeval cputime-cpt];
        logstruct.condestDF = condest(DF);
    end

    cpt = cputime;
    F=OCMATCONT.operatoreq(tmesh,coeff,tangent);  % new residual
    if log 
        logstruct.runtime.Feval = [logstruct.runtime.Feval cputime-cpt]; 
    end
    fcount = fcount +1;            % increase number of function evaluations
end

cpt=cputime;
[L,U]=lu(DF);                  % LU-decomposition of DF
if log 
    logstruct.runtime.lu = [logstruct.runtime.lu cputime-cpt]; 
end

lastwarn('');                  % initialize warning state

cpt=cputime;
if OCMATCONT.OPTIONS.continuationmethod==1 % Moore Penrose
    % find the Newton direction
    D=U\(L\([-F R]));         % Newton correction
    delta_coeff=D(:,1);
    tangent=D(:,2);
    tangent=tangent/norm(tangent);
else % Arclength
    delta_coeff = U\(L\(-F ));
end
if log
    logstruct.runtime.resubstitute = [logstruct.runtime.resubstitute cputime-cpt];
end

if length(lastwarn)            % Exit if DF is singular
    %   error(' System matrix is close to singular. Try a refined mesh or another initial approximation.');
end

G0=norm(delta_coeff);            % Norm of the Newton correction

new_coeff=coeff+delta_coeff;     % new approximation

% ***** In case we have very good initial approximation (e.g. from a previous mesh),
% ***** tolerances might be satisfied even now and we reduce computational
% ***** costs to the absolute minimum of one F- and one DF-evaluation
TolFactor = check_tolerances(new_coeff,delta_coeff);
if TolFactor < 1
    DOC = [];  % dummy outputs
    G = [];
    simplified_delta_coeff=[];
    simplified_tangent=[];
    if log
        logstruct.tol_factor=[logstruct.tol_factor TolFactor];
        logstruct.G0= [logstruct.G0 G0];
    end
    DF=OCMATCONT.frechetder(tmesh,new_coeff,tangent);
    DF=[DF;tangent(:).'];
    [L,U] = lu(DF);
    tangent=U\(L\(R));
    tangent=tangent/norm(tangent);
    return
end

% ********** Check Monotonicity Condition for lambda = 1 and lambda = lambda_min
cpt = cputime;
F=OCMATCONT.operatoreq(tmesh,new_coeff,tangent);  % new residual
if log 
    logstruct.runtime.Feval = [logstruct.runtime.Feval cputime-cpt]; 
end
fcount = fcount +1;             % increase number of function evaluations

cpt = cputime;
if OCMATCONT.OPTIONS.continuationmethod==1 % Moore Penrose
    % find the Newton direction
    D=U\(L\([-F R]));         % Newton correction
    simplified_delta_coeff=D(:,1);
    simplified_tangent=D(:,2);
    simplified_tangent=simplified_tangent/norm(simplified_tangent);
else % Arclength
    simplified_delta_coeff = U\(L\(-F));            % simplified delta_coeff
end
if log 
    logstruct.runtime.resubstitute = [logstruct.runtime.resubstitute cputime-cpt]; 
end

G = norm(simplified_delta_coeff);

if G <(1-OCMATCONT.OPTIONS.lambdamin /2) * G0 % Approximation improved
    if G<OCMATCONT.OPTIONS.switchtoffnfactor*G0         % Approximation improved significantly
        DOC = 1;             % -> Fast Frozen Newton

        % See if tolerances are satisfied
        TolFactor = check_tolerances(new_coeff,simplified_delta_coeff);
        if TolFactor < 1
            % Perform one last step towards the solution
            new_coeff= new_coeff+simplified_delta_coeff;

            if log
                logstruct.tol_factor = [logstruct.tol_factor TolFactor];
                logstruct.G0 = [logstruct.G0 G0];
                logstruct.G = [logstruct.G G];
            end
            DF=OCMATCONT.frechetder(tmesh,new_coeff,tangent);
            DF=[DF;simplified_tangent(:).'];
            [L,U] = lu(DF);
            tangent=U\(L\(R));
            tangent=tangent/norm(tangent);
            return
        end
    else % Approximation improved slightly
        DOC = 2;             % -> Predictor Corrector Line Search
        if log
            logstruct.jac_update  = [logstruct.jac_update  0];   % Use Jacobian computed here
        end
    end
else % try another Newton step with smallest possible lambda

    new_coeff=coeff+OCMATCONT.OPTIONS.lambdamin*delta_coeff;      % new approximation
    cpt = cputime;
    Fmin=OCMATCONT.operatoreq(tmesh,new_coeff,tangent);   % new residual
    if log 
        logstruct.runtime.Feval = [logstruct.runtime.Feval cputime-cpt];
    end
    fcount = fcount +1;
    cpt = cputime;
    simplified_delta_coeff_min = U\(L\(-Fmin));
    if log
        logstruct.runtime.resubstitute = [logstruct.runtime.resubstitute cputime-cpt];
    end

    Gmin = norm(simplified_delta_coeff_min);

    if Gmin < (1-OCMATCONT.OPTIONS.lambdamin /2) * G0 % Approximation improved
        DOC = 2; % Predictor Corrector Line Search
        if log
            logstruct.jac_update  = [logstruct.jac_update  0];   % Use Jacobian computed here
        end
    else
        DOC = 3; % -> Trust Region Method
    end
end


%********************************************************************************
%******************************* CHECK TOLERANCES *******************************
%********************************************************************************

function TolFactor = check_tolerances(coeff,delta)
% Calculates the factor by which delta has to be reduced to satisfy the
% tolerances.
% out <= 1   => Tolerances satisfied
% out >  1   => Tolerances not satisfied
global OCMATCONT


AbsTol=OCMATCONT.OPTIONS.newtonabstol;
RelTol=OCMATCONT.OPTIONS.newtonreltol;

% ***** Determine Dimension of solution array coeff
[d p1 N] = size(coeff);

% ***** Extract solution components at mesh points (last meshpoint (b) not included
InfNorm_y = max(abs(coeff));

% ***** Extract last correction at mesh points
InfNorm_delta = max(abs(delta));

% ***** Componentwise check if tolerances are satisfied
TolFactor = max(InfNorm_delta./(AbsTol+InfNorm_y*RelTol));
