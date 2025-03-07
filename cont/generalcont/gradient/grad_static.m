function [extremal,graditer,norm_dx,R]=grad_static(extremal,par)
global OCGRADCONT OCGRADSOL

% cfv ... constraint function value
x = extremal.y;
R=OCGRADCONT.OPTIONS.R0; % penalty factor for violation


norm_dx=OCGRADCONT.OPTIONS.gradtol+1;
graditer = 1; % design point number.  graditer = 0 is initial design point.

% enter into a loop if initial design does not satisfy the condition of
% length of the gradient being smaller than epsilon.
jmax=-inf;
while 1

    % Values of the objective function, constraints, and gradients at the
    % current design point.
    Vnew=OCGRADSOL.objective(x,par);

    h_value=OCGRADSOL.equalityconstraint(x,par);
    if OCGRADSOL.equalitynum
        h_value_max=abs(min(0,min(h_value)));
        grad_h_b_value=OCGRADSOL.gradientec(x,par);
    else
        h_value_max=0;
        grad_h_b_value=[];
    end

    g_value=OCGRADSOL.inequalityconstraint(x,par);
    if OCGRADSOL.inequalitynum
        g_value_max=abs(min(0,min(g_value)));
        grad_g_b_value=OCGRADSOL.gradientic(x,par);
    else
        g_value_max=0;
        grad_g_b_value=[];
    end
    cfv=max(0,max(h_value_max,g_value_max)); % maximum violation of constraint


    % solve QP subproblem
    switch OCGRADCONT.OPTIONS.gradientmappingmethod
        case 'explicit'
            [dx,lmec_qp,lmiec_qp]=optqpsol(x,par); % dx ... admissible descent direction
        case 'numeric'
            c=-OCGRADSOL.gradientobjective(x,par);
            b=[h_value;g_value];
            A=[grad_h_b_value;grad_g_b_value];
            [dx,fval,exitflag,output,lm]=quadprog(diag(OCGRADSOL.d),c,A,b,[],[],[],[],zeros(OCGRADSOL.statenum,1),OCGRADCONT.OPTIONS.optimset);
            lmiec_qp=lm.ineqlin;
            lmec_qp=lm.eqlin;
            %
            %             [dx0,lmec_qp0,lmiec_qp0]=optqpsol(x,par); % dx ... admissible descent direction
            %             if abs(dx-dx0)>1e-6 | length(dx0)>1
            %                 dx0
            %             end
        case 'nesterov'
            [dx,lmec_qp,lmiec_qp]=nesterovsol(x,par);
    end
    norm_dx = norm(dx);

    if norm_dx<OCGRADCONT.OPTIONS.gradtol % if admissible change in admissible direction is small stop
        break
    end

    r = 0;
    if OCGRADSOL.equalitynum > 0
        r = r + sum(abs(lmec_qp));
    end
    if OCGRADSOL.inequalitynum > 0
        r = r + sum(lmiec_qp);
    end

    R = max(R, r);


    % Inexact line search
    step=OCGRADCONT.OPTIONS.initlinestepwidth;
    phi = Vnew-R*cfv;
    beta = OCGRADCONT.OPTIONS.gamma*norm_dx.^2;

    linesearchiter = 0;
    cond = false;
    while cond == false && linesearchiter < OCGRADCONT.OPTIONS.maxlinesearchiter

        t = (step).^(linesearchiter);
        %
        x_next=x+t*dx;
        Vnew_next =  OCGRADSOL.objective(x_next,par);

        h_value=OCGRADSOL.equalityconstraint(x_next,par);
        if OCGRADSOL.equalitynum
            h_value_max=abs(min(0,min(h_value)));
        else
            h_value_max=0;
        end

        g_value=OCGRADSOL.inequalityconstraint(x_next,par);
        if OCGRADSOL.inequalitynum
            g_value_max=abs(min(0,min(g_value)));
        else
            g_value_max=0;
        end

        cfv_next = max(0,max(h_value_max,g_value_max));

        phi_next=Vnew_next-R*cfv_next;

        if phi_next>phi+t*beta
            cond = true;
            x = x_next;
        end

        linesearchiter=linesearchiter+1;
    end
    jmax=max(jmax,linesearchiter);
    graditer = graditer + 1;
    if graditer>OCGRADCONT.OPTIONS.maxgraditer
        break
    end
end
extremal.y=x;
%extremal.tangent=[];
extremal.objectivevalue=Vnew;


function [dx,lmec,lmiec]=optqpsol(x,par)
global OCGRADCONT OCGRADSOL
stateidx=1:OCGRADSOL.statenum;
if OCGRADSOL.equalitynum
    lmecidx=(1:OCGRADSOL.equalitynum)+stateidx(end);
else
    lmecidx=stateidx(end);
end
if OCGRADSOL.inequalitynum
    lmiecidx=(1:OCGRADSOL.inequalitynum)+lmecidx(end);
    slvidx=(1:OCGRADSOL.inequalitynum)+lmiecidx(end);
else
    lmiecidx=0;
    slvidx=0;
end

kkp=OCGRADSOL.qpsolution(x,par);

kkp(:,sum(abs(imag(kkp)))>0)=[];
%kkp(:,OCGRADSOL.inequalityconstraint(x,par)-kkp(1,:)>0)=[];
dx=kkp(stateidx,:);
[dx,idx]=uniqueoc(dx);

if OCGRADSOL.equalitynum
    lmec=kkp(lmecidx,idx);
else
    lmec=[];
end
if OCGRADSOL.inequalitynum
    lmiec=kkp(lmiecidx,idx);
    slv=kkp(slvidx,idx);
else
    lmiec=[];
    slv=[];
end
idx=find(sum(abs(lmiec<0))>0);
dx(:,idx)=[];
if OCGRADSOL.equalitynum
    lmec(:,idx)=[];
end
if OCGRADSOL.inequalitynum
    lmiec(:,idx)=[];
end


function [dx,lmec,lmiec]=nesterovsol(x,par)
global OCGRADCONT OCGRADSOL
stateidx=1:OCGRADSOL.statenum;
if OCGRADSOL.equalitynum
    lmecidx=(1:OCGRADSOL.equalitynum)+stateidx(end);
else
    lmecidx=stateidx(end);
end
if OCGRADSOL.inequalitynum
    lmiecidx=(1:OCGRADSOL.inequalitynum)+lmecidx(end);
    slvidx=(1:OCGRADSOL.inequalitynum)+lmiecidx(end);
else
    lmiecidx=0;
    slvidx=0;
end

kkp=OCGRADSOL.nesterovsolution(x,par,OCGRADCONT.OPTIONS.gradientmappinggamma);

kkp(:,sum(abs(imag(kkp)))>0)=[];

xq=kkp(stateidx,:);
[xq,idx]=uniqueoc(xq);

if OCGRADSOL.equalitynum
    lmec=kkp(lmecidx,idx);
else
    lmec=[];
end
if OCGRADSOL.inequalitynum
    lmiec=kkp(lmiecidx,idx);
    slv=kkp(slvidx,idx);
else
    lmiec=[];
    slv=[];
end
idx=find(sum(abs(lmiec<0))>0);
xq(:,idx)=[];
if OCGRADSOL.equalitynum
    lmec(:,idx)=[];
end
if OCGRADSOL.inequalitynum
    lmiec(:,idx)=[];
end

dx=OCGRADCONT.OPTIONS.gradientmappinggamma*(xq-x);