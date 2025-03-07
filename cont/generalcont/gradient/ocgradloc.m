function [extremal,graditer,max_lineiter,statinfo]=ocgradloc(extremal,par)
% solves a (local) optimal control problem using the gradient method
%
% Y ... vector consisting of (local) states and objective value
% t ... time grid
% par ... modelparameter
% v ... (local) control vector
% statedynamicsext: statedynamics+objectivefunction

global OCGRADCONT OCGRADSOL

t=extremal.t;
T=extremal.timehorizon; % time horizon
Y0=[extremal.y(:,1);0];
cst_YT=extremal.cst_y(:,end);

v=extremal.v;
graditer=0;
linestepwidth=OCGRADCONT.OPTIONS.initlinestepwidth; % stepwidth for the linesearch
statinfo=struct('lineiter',{},'linestepwidth',{});

if OCGRADCONT.cconstraint
    if OCGRADCONT.mcconstraint
        adm=0;
        ctr=0;
        while ~adm
            ctr=ctr+1;
            Y=OCGRADCONT.gradientodesolver(1,OCGRADSOL.statedynamicsext,t,Y0,v,par);
            adm=OCGRADSOL.admissible(t,Y,v,par);
            if adm || ctr>5
                break
            end
            v=OCGRADSOL.projectionlocal(Y,v,par);
        end
    else
        if OCGRADCONT.OPTIONS.projection
            v=OCGRADSOL.projectionlocal(v,par);
        end
        Y=OCGRADCONT.gradientodesolver(1,OCGRADSOL.statedynamicsext,t,Y0,v,par);
    end
else
    Y=OCGRADCONT.gradientodesolver(1,OCGRADSOL.statedynamicsext,t,Y0,v,par);
end
max_lineiter=-inf;
while 1
    graditer=graditer+1;
    %oval=objectivevalue(T,Y,par);

    cst_y0=endconstraint(T,Y,par,cst_YT);
    cst_y=OCGRADCONT.gradientodesolver(-1,OCGRADSOL.costatedynamics,t,cst_y0,[Y(OCGRADCONT.statecoordinate.concentrated,:);v],par);

    % calculate the new search direction
    dv=gradientcontrol(t,Y(OCGRADCONT.statecoordinate.concentrated,:),cst_y,v,par);

    v0=v;
    switch OCGRADCONT.linesearchmethod
        case 'min'
            [Y,v,linestepwidth,impr,lineiter]=linesearchloc_min(t,Y,v,dv,par,linestepwidth);
        case 'armijo'
            [Y,v,linestepwidth0,impr,lineiter]=linesearchloc_armijo(t,Y,v,dv,par,linestepwidth);
    end
    if OCGRADCONT.OPTIONS.globalcorrector
        v=OCGRADSOL.globalcorrector(t,Y,v,cst_y,par);
    end
    dv0=v0-v;
    if OCGRADCONT.control_num.concentrated>1
        %norm_dv=max(sqrt(sum(dv.*dv)));
        norm_dv=sqrt(sum(sum(dv0.*dv0)));
    else
        %norm_dv=max(sqrt(dv.*dv));
        norm_dv=sqrt(sum(dv0.*dv0));
    end

    if norm_dv<OCGRADCONT.OPTIONS.gradtol || ~impr || graditer>OCGRADCONT.OPTIONS.maxgraditer || linestepwidth<=OCGRADCONT.OPTIONS.minlinestepwidth% if admissible change in admissible direction is small stop
        break
    end
    max_lineiter=max(max_lineiter,lineiter);
    if OCGRADCONT.OPTIONS.printgradientinfo
        statinfo(graditer).lineiter=lineiter;
        statinfo(graditer).linestepwidth=linestepwidth;
        statinfo(graditer).norm_dv=norm_dv;
    end
end


% if linestepwidth<=OCGRADCONT.OPTIONS.minlinestepwidth
%     ocmatmsg('Gradient algorithm may not converge, minimum linestep width reached.');
% end
Y=OCGRADCONT.gradientodesolver(1,OCGRADSOL.statedynamicsext,t,Y0,v,par);
cst_y0=endconstraint(T,Y,par,cst_YT);
cst_y=OCGRADCONT.gradientodesolver(-1,OCGRADSOL.costatedynamics,t,cst_y0,[Y(OCGRADCONT.statecoordinate.concentrated,:);v],par);
dv=gradientcontrol(t,Y(OCGRADCONT.statecoordinate.concentrated,:),cst_y,v,par);

extremal.v=v;
extremal.dv=dv;
extremal.y=Y(OCGRADCONT.statecoordinate.concentrated,:);
extremal.o=Y(OCGRADCONT.objectivecoordinate,:);
extremal.objectivevalue=objectivevalue(T,Y,par,v);
extremal.cst_y=cst_y;
extremal.coeff=[v(:);extremal.y(:,1);extremal.cst_y(:,OCGRADCONT.TIMEMESH.num)];

function dv=gradientcontrol(t,y,cst_y,v,par)
global OCGRADCONT OCGRADSOL

dv=OCGRADSOL.gradhamiltonian(t,y,v,cst_y,par);
switch OCGRADCONT.OPTIONS.gradientmappingmethod
    case 'nesterov'
        dv=nesterovmapping(t,y,cst_y,v,par,dv);

%                 [xq0,lm0]=OCGRADSOL.explicitgradientcontrol(t,y,v,dv,par,OCGRADCONT.OPTIONS.gradientmappinggamma);
%                 dv0=OCGRADCONT.OPTIONS.gradientmappinggamma*(xq0-v);
    case 'explicit'
        try
            xq=OCGRADSOL.explicitgradientcontrol(t,y,v,dv,par,OCGRADCONT.OPTIONS.gradientmappinggamma);
            dv=OCGRADCONT.OPTIONS.gradientmappinggamma*(xq-v);
        catch
            dv=nesterovmapping(t,y,cst_y,v,par,dv);
        end
        %         % for testing
        %         [lb,ub]=OCGRADSOL.controlbounds(t,y,par);
        %         A=[];
        %         b=[];
        %
        %         flin=-dv-OCGRADCONT.OPTIONS.gamma*v;
%         flin=flin(:);
%         fquad=OCGRADCONT.Nesterov.QuadraticMatrix;
%         xqtst=quadprog(fquad,flin,A,b,[],[],lb,ub,-dv,OCGRADCONT.OPTIONS.optimset).';
%         dvtst=OCGRADCONT.OPTIONS.gamma*(xqtst-v);
        
end

function dv=nesterovmapping(t,y,cst_y,v,par,dv)
global OCGRADCONT OCGRADSOL

[lb,ub,A,b]=OCGRADSOL.controlbounds(t,y,par);
flin=-dv-OCGRADCONT.OPTIONS.gradientmappinggamma*v;
flin=flin(:);
fquad=OCGRADCONT.Nesterov.QuadraticMatrix;
xq=quadprog(fquad,flin,A,b,[],[],lb,ub,[],OCGRADCONT.OPTIONS.optimset);
xq=reshape(xq,[],length(t));
dv=OCGRADCONT.OPTIONS.gradientmappinggamma*(xq-v);

function [Y_best,v_best,step_best,impr,lineiter]=linesearchloc_min(t,Y,v,dv,par,step)
global OCGRADCONT OCGRADSOL

% Y consists of states, costates and objective path
impr=0;

T=t(OCGRADCONT.TIMEMESH.num);
v0=v;
step0=step;
%oval=Y(end,end);
lineiter=0;
while impr==0 && lineiter<=OCGRADCONT.OPTIONS.maxlinesearchiter && step>OCGRADCONT.OPTIONS.minlinestepwidth
    step = min(step,OCGRADCONT.OPTIONS.maxlinestepwidth);
    y0=objectivevalue(T,Y,par,v0);
    v=v0+step*dv;
    if OCGRADCONT.cconstraint
        if OCGRADCONT.mcconstraint
            adm=0;
            ctr=0;
            while ~adm
                ctr=ctr+1;
                Y_new=OCGRADCONT.gradientodesolver(1,OCGRADSOL.statedynamicsext,t,Y(:,1),v,par);
                adm=OCGRADSOL.admissible(t,Y_new,v,par);
                if adm || ctr>5
                    break
                end
                v=OCGRADSOL.projectionlocal(Y_new,v,par);
            end
        else
            if OCGRADCONT.OPTIONS.projection
                v=OCGRADSOL.projectionlocal(v,par);
            end
            Y_new=OCGRADCONT.gradientodesolver(1,OCGRADSOL.statedynamicsext,t,Y(:,1),v,par);
        end
    else
        Y_new=OCGRADCONT.gradientodesolver(1,OCGRADSOL.statedynamicsext,t,Y(:,1),v,par);
    end
    y1=objectivevalue(T,Y_new,par,v);

    lineiter=lineiter+1;
    x1=step;

    if (y1>y0) %%% 1.1
        impr = 1;
        v_best=v;
        Y_best=Y_new;
        step_best=step;

        step = (1+OCGRADCONT.OPTIONS.linestep_increase)*step;
        x2 = step;
        v=v0+step*dv;
        if OCGRADCONT.cconstraint
            if OCGRADCONT.mcconstraint
                adm=0;
                ctr=0;
                while ~adm
                    ctr=ctr+1;
                    Y_new=OCGRADCONT.gradientodesolver(1,OCGRADSOL.statedynamicsext,t,Y(:,1),v,par);
                    adm=OCGRADSOL.admissible(t,Y_new,v,par);
                    if adm || ctr>5
                        break
                    end
                    v=OCGRADSOL.projectionlocal(Y_new,v,par);
                end
            else
                if OCGRADCONT.OPTIONS.projection
                    v=OCGRADSOL.projectionlocal(v,par);
                end
                Y_new=OCGRADCONT.gradientodesolver(1,OCGRADSOL.statedynamicsext,t,Y(:,1),v,par);
            end
        else
            Y_new=OCGRADCONT.gradientodesolver(1,OCGRADSOL.statedynamicsext,t,Y(:,1),v,par);
        end
        oval_new=objectivevalue(T,Y_new,par,v);
        lineiter=lineiter+1;

        y2=oval_new;

        if (y2>y1) %%% 1.1.1
            v_best=v;
            Y_best=Y_new;
            step_best=step;
            if (x1*y2+x2*y0-(x1+x2)*y1<=0) %%% 1.1.1.1      %Check if curve is convex
                % Curve is convex --> second solution is the best
            else %%% 1.1.1.1
                % Curve is concave --> stepsize for maximum of quad. approximation
                step=((x2^2-x1^2)*y0+y2*x1^2-y1*x2^2)/(2*(y2*x1-y1*x2+y0*(x2-x1)));
                step = min(step,OCGRADCONT.OPTIONS.maxlinestepwidth);

                % New controls + new states + new objective function
                v=v0+step*dv;
                if OCGRADCONT.cconstraint
                    if OCGRADCONT.mcconstraint
                        adm=0;
                        ctr=0;
                        while ~adm
                            ctr=ctr+1;
                            Y_new=OCGRADCONT.gradientodesolver(1,OCGRADSOL.statedynamicsext,t,Y(:,1),v,par);
                            adm=OCGRADSOL.admissible(t,Y_new,v,par);
                            if adm || ctr>5
                                break
                            end
                            v=OCGRADSOL.projectionlocal(Y_new,v,par);
                        end
                    else
                        if OCGRADCONT.OPTIONS.projection
                            v=OCGRADSOL.projectionlocal(v,par);
                        end
                        Y_new=OCGRADCONT.gradientodesolver(1,OCGRADSOL.statedynamicsext,t,Y(:,1),v,par);
                    end
                else
                    Y_new=OCGRADCONT.gradientodesolver(1,OCGRADSOL.statedynamicsext,t,Y(:,1),v,par);
                end
                oval_new=objectivevalue(T,Y_new,par,v);
                lineiter=lineiter+1;

                if oval_new>y2 %%% 1.1.1.1.1
                    % Save currently best values
                    v_best=v;
                    Y_best=Y_new;
                    step_best=step;

                end %%% 1.1.1.1.1
            end %%% 1.1.1.1
        else %%% 1.1.1         %(y2 <= y1)
            step = ((x2^2-x1^2)*y0+y2*x1^2-y1*x2^2)/(2*(y2*x1-y1*x2+y0*(x2-x1)));
            step = min(step,OCGRADCONT.OPTIONS.maxlinestepwidth);

            % New controls + new states + new objective function
            v=v0+step*dv;
            if OCGRADCONT.cconstraint
                if OCGRADCONT.mcconstraint
                    adm=0;
                    ctr=0;
                    while ~adm
                        ctr=ctr+1;
                        Y_new=OCGRADCONT.gradientodesolver(1,OCGRADSOL.statedynamicsext,t,Y(:,1),v,par);
                        adm=OCGRADSOL.admissible(t,Y_new,v,par);
                        if adm || ctr>5
                            break
                        end
                        v=OCGRADSOL.projectionlocal(Y_new,v,par);
                    end
                else
                    if OCGRADCONT.OPTIONS.projection
                        v=OCGRADSOL.projectionlocal(v,par);
                    end
                    Y_new=OCGRADCONT.gradientodesolver(1,OCGRADSOL.statedynamicsext,t,Y(:,1),v,par);
                end
            else
                Y_new=OCGRADCONT.gradientodesolver(1,OCGRADSOL.statedynamicsext,t,Y(:,1),v,par);
            end
            oval_new=objectivevalue(T,Y_new,par,v);
            lineiter=lineiter+1;

            if (oval_new>y1) %%% 1.1.1.2

                % Save currently best values
                v_best=v;
                Y_best=Y_new;
                step_best=step;
            end %%% 1.1.1.2
        end %%% 1.1.1
    else %%% 1.1          (y1 <= y0)
        step = 0.5*step;
        x2 = step;
        v=v0+step*dv;
        if OCGRADCONT.cconstraint
            if OCGRADCONT.mcconstraint
                adm=0;
                ctr=0;
                while ~adm
                    ctr=ctr+1;
                    Y_new=OCGRADCONT.gradientodesolver(1,OCGRADSOL.statedynamicsext,t,Y(:,1),v,par);
                    adm=OCGRADSOL.admissible(t,Y_new,v,par);
                    if adm || ctr>5
                        break
                    end
                    v=OCGRADSOL.projectionlocal(Y_new,v,par);
                end
            else
                if OCGRADCONT.OPTIONS.projection
                    v=OCGRADSOL.projectionlocal(v,par);
                end
                Y_new=OCGRADCONT.gradientodesolver(1,OCGRADSOL.statedynamicsext,t,Y(:,1),v,par);
            end
        else
            Y_new=OCGRADCONT.gradientodesolver(1,OCGRADSOL.statedynamicsext,t,Y(:,1),v,par);
        end
        oval_new=objectivevalue(T,Y_new,par,v);
        lineiter=lineiter+1;

        y2=oval_new;
        if (y2>y0) %%% 1.1.2

            impr = 1;

            % Save currently best values
            v_best=v;
            Y_best=Y_new;
            step_best=step;

            % Curve is concave --> stepsize for maximum of quad. approximation
            step=((x2^2-x1^2)*y0+y2*x1^2-y1*x2^2)/(2*(y2*x1-y1*x2+y0*(x2-x1)));
            step = min(step,OCGRADCONT.OPTIONS.maxlinestepwidth);

            % New controls + new states + new objective function
            v=v0+step*dv;
            if OCGRADCONT.cconstraint
                if OCGRADCONT.mcconstraint
                    adm=0;
                    ctr=0;
                    while ~adm
                        ctr=ctr+1;
                        Y_new=OCGRADCONT.gradientodesolver(1,OCGRADSOL.statedynamicsext,t,Y(:,1),v,par);
                        adm=OCGRADSOL.admissible(t,Y_new,v,par);
                        if adm || ctr>5
                            break
                        end
                        v=OCGRADSOL.projectionlocal(Y_new,v,par);
                    end
                else
                    if OCGRADCONT.OPTIONS.projection
                        v=OCGRADSOL.projectionlocal(v,par);
                    end
                    Y_new=OCGRADCONT.gradientodesolver(1,OCGRADSOL.statedynamicsext,t,Y(:,1),v,par);
                end
            else
                Y_new=OCGRADCONT.gradientodesolver(1,OCGRADSOL.statedynamicsext,t,Y(:,1),v,par);
            end
            oval_new=objectivevalue(T,Y_new,par,v);
            lineiter=lineiter+1;

            if (oval_new>y2) %%% 1.1.2.1

                % Save currently best values
                v_best=v;
                Y_best=Y_new;
                step_best=step;
            end %%% 1.1.2.1
        elseif (y2>0.5*(y0+y1)) % && y2 < (3*y0+y1)/4   %%% 1.1.2
            step=((x2^2-x1^2)*y0+y2*x1^2-y1*x2^2)/(2*(y2*x1-y1*x2+y0*(x2-x1)));
            step = min(step,OCGRADCONT.OPTIONS.maxlinestepwidth);
            if (step > 0) %%% 1.1.2.2

                % New controls + new states + new objective function
                v=v0+step*dv;
                if OCGRADCONT.cconstraint
                    if OCGRADCONT.mcconstraint
                        adm=0;
                        ctr=0;
                        while ~adm
                            ctr=ctr+1;
                            Y_new=OCGRADCONT.gradientodesolver(1,OCGRADSOL.statedynamicsext,t,Y(:,1),v,par);
                            adm=OCGRADSOL.admissible(t,Y_new,v,par);
                            if adm || ctr>5
                                break
                            end
                            v=OCGRADSOL.projectionlocal(Y_new,v,par);
                        end
                    else
                        if OCGRADCONT.OPTIONS.projection
                            v=OCGRADSOL.projectionlocal(v,par);
                        end
                        Y_new=OCGRADCONT.gradientodesolver(1,OCGRADSOL.statedynamicsext,t,Y(:,1),v,par);
                    end
                else
                    Y_new=OCGRADCONT.gradientodesolver(1,OCGRADSOL.statedynamicsext,t,Y(:,1),v,par);
                end
                oval_new=objectivevalue(T,Y_new,par,v);
                lineiter=lineiter+1;

                if (oval_new>=y0) %%% 1.1.2.2.1
                    impr = 1;

                    % Save currently best values
                    v_best=v;
                    Y_best=Y_new;
                    step_best=step;

                end %%% 1.1.2.2.1
            end %%% 1.1.2.2
        end %%% 1.1.2
    end
end
if impr == 0
    % If there is no improvement, return original values
    v_best=v0;
    Y_best=Y;
    step_best=step0;
    % only step-size gets adjusted
    % --> new try from the beginning with the new step-size
end

function [Y_best,v_best,step_best,impr,lineiter]=linesearchloc_armijo(t,Y,v0,dv,par,step)
global OCGRADCONT OCGRADSOL

% Y consists of states, costates and objective path
impr=false;
T=t(OCGRADCONT.TIMEMESH.num);
beta=OCGRADCONT.OPTIONS.gamma;
mu=OCGRADCONT.OPTIONS.mu;
oval0=objectivevalue(T,Y,par,v0);

step0=step;

lineiter=1;
while lineiter<=OCGRADCONT.OPTIONS.maxlinesearchiter

    v_new=v0+step*dv;
    
    if OCGRADCONT.cconstraint
        if OCGRADCONT.mcconstraint
            adm=0;
            ctr=0;
            while ~adm
                ctr=ctr+1;
                Y_new=OCGRADCONT.gradientodesolver(1,OCGRADSOL.statedynamicsext,t,Y(:,1),v_new,par);
                adm=OCGRADSOL.admissible(t,Y_new,v_new,par);
                if adm || ctr>5
                    break
                end
                v_new=OCGRADSOL.projectionlocal(Y_new,v_new,par);
            end
        else
            if OCGRADCONT.OPTIONS.projection
                v_new=OCGRADSOL.projectionlocal(v_new,par);
            end
            Y_new=OCGRADCONT.gradientodesolver(1,OCGRADSOL.statedynamicsext,t,Y(:,1),v_new,par);
        end
    else
        Y_new=OCGRADCONT.gradientodesolver(1,OCGRADSOL.statedynamicsext,t,Y(:,1),v_new,par);
    end
    oval_new=objectivevalue(T,Y_new,par,v_new);

    if oval_new>=oval0+mu*sum(sum(dv.*(v_new-v0)))
        impr=true;
        v_best=v_new;
        Y_best=Y_new;
        step_best=step;
        break
    end
    lineiter=lineiter+1;
    step=step*beta;
end
if ~impr
    % If there is no improvement, return original values
    v_best=v0;
    Y_best=Y;
    step_best=step0;
    % only step-size gets adjusted
    % --> new try from the beginning with the new step-size
end

function oval=objectivevalue(T,Y,par,v)
global OCGRADCONT OCGRADSOL

if OCGRADCONT.infinitetimeendconditions
    oval=Y(OCGRADCONT.objectivecoordinate,OCGRADCONT.TIMEMESH.num)+OCGRADSOL.asymptoticobjectivefunction(T,Y(OCGRADCONT.statecoordinate.concentrated,OCGRADCONT.TIMEMESH.num),v(:,OCGRADCONT.TIMEMESH.num),par);
else
    oval=Y(OCGRADCONT.objectivecoordinate,OCGRADCONT.TIMEMESH.num)+OCGRADSOL.salvagevalue(T,Y(OCGRADCONT.statecoordinate.concentrated,OCGRADCONT.TIMEMESH.num),par);
end
function cst_y=endconstraint(T,Y,par,cst_YT)
global OCGRADCONT OCGRADSOL

if OCGRADCONT.infinitetimeendconditions
    cst_y=OCGRADSOL.asymptotictransversaltycondition(T,Y(OCGRADCONT.statecoordinate.concentrated,OCGRADCONT.TIMEMESH.num),cst_YT,par);
else
    cst_y=OCGRADSOL.transversaltycondition(T,Y(OCGRADCONT.statecoordinate.concentrated,OCGRADCONT.TIMEMESH.num),cst_YT,par);
end

