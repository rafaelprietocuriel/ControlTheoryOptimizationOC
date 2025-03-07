function [extremal,graditer,max_lineiter]=ocgradlocp(extremal,par)
% solves a (local) optimal control problem using the gradient method
%
% Y ... vector consisting of (local) states and objective value
% t ... time grid
% par ... modelparameter
% v ... (local) control vector
% statedynamicsext: statedynamics+objectivefunction

global OCGRADCONT OCGRADSOL

%norm_dv=OCGRADCONT.OPTIONS.gradtol+1;
R=OCGRADCONT.OPTIONS.R0; % penalty factor for violation

t=extremal.t;
T=extremal.timehorizon; % time horizon
Y0=[extremal.y(:,1);0];
cst_YT=extremal.cst_y(:,end);
v=extremal.v;

Y=OCGRADCONT.gradientodesolver(1,OCGRADSOL.statedynamicsext,t,Y0,v,par);

graditer=0;
impr=1; % flag for improvement of the objective value, repeat gradient steps until no improvement is possible

max_lineiter=-inf;
while 1
    graditer=graditer+1;

    if OCGRADCONT.cconstraint
        cstr=OCGRADSOL.constraint(t,Y,v,par);
        cstr(cstr>0)=0;
    else
        cstr=0;
    end

    cfv=max(0,abs(min(cstr(:)))); % maximum violation of constraint

    cst_y0=endconstraint(T,Y,par,cst_YT);
    cst_y=OCGRADCONT.gradientodesolver(-1,OCGRADSOL.costatedynamics,t,cst_y0,[Y(OCGRADCONT.statecoordinate.concentrated,:);v],par);

    % calculate the new search direction
    [dv,lm]=gradientcontrol(t,Y(OCGRADCONT.statecoordinate.concentrated,:),cst_y,v,par);

    if OCGRADCONT.control_num.concentrated>1
        norm_dv=max(sqrt(sum(dv.*dv)));
    else
        norm_dv=max(sqrt(dv.*dv));
    end

    if norm_dv<OCGRADCONT.OPTIONS.gradtol || ~impr || graditer>OCGRADCONT.OPTIONS.maxgraditer% if admissible change in admissible direction is small stop
        break
    end
    
    %
    if ~isempty(lm)
        lm(lm>0)=0;
        R=max(R,max(abs(min(lm))));
    end
    
    [Y,v,impr,lineiter]=linesearchloc(t,Y,v,dv,par,R,cfv,norm_dv);
    if OCGRADCONT.OPTIONS.globalcorrector
        v=OCGRADSOL.globalcorrector(t,Y,v,cst_y,par);
    end
    max_lineiter=max(max_lineiter,lineiter);
end

if graditer>OCGRADCONT.OPTIONS.maxgraditer
    ocmatmsg('Gradient algorithm may not converge, number of maximum gradient steps exceeded.');
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
%extremal.oinf=OCGRADSOL.asymptoticobjectivevalue(T,Y(OCGRADCONT.statecoordinate.concentrated,1),v(:,1),cst_y(:,1),par);

%extremal.improvement=impr; % flag that expresses if this solution is an improvement to the original solution

function [dv,lm]=gradientcontrol(t,y,cst_y,v,par)
global OCGRADCONT OCGRADSOL

dv=OCGRADSOL.gradhamiltonian(t,y,v,cst_y,par);
switch OCGRADCONT.OPTIONS.gradientmappingmethod
    case 'nesterov'
        [xq0,lm0]=OCGRADSOL.explicitgradientcontrol(t,y,v,dv,par,OCGRADCONT.OPTIONS.gradientmappinggamma);
        dv0=OCGRADCONT.OPTIONS.gradientmappinggamma*(xq0-v);

        [dv,lm]=nesterovmapping(t,y,cst_y,v,par,dv);
    case 'explicit'
        try
            [xq,lm]=OCGRADSOL.explicitgradientcontrol(t,y,v,dv,par,OCGRADCONT.OPTIONS.gradientmappinggamma);
            dv=OCGRADCONT.OPTIONS.gradientmappinggamma*(xq-v);
        catch
            [dv,lm]=nesterovmapping(t,y,cst_y,v,par,dv);
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
        
    otherwise
        lm=[];
end

function [dv,lm]=nesterovmapping(t,y,cst_y,v,par,dv)
global OCGRADCONT OCGRADSOL

[lb,ub,A,b]=OCGRADSOL.controlbounds(t,y,par);
flin=-dv-OCGRADCONT.OPTIONS.gradientmappinggamma*v;
flin=flin(:);
fquad=OCGRADCONT.Nesterov.QuadraticMatrix;
[xq,fval,exitflag,output,tmplm]=quadprog(fquad,flin,A,b,[],[],lb,ub,[],OCGRADCONT.OPTIONS.optimset);
xq=reshape(xq,[],length(t));
if ~isempty(lb)
    tmplm.lower=reshape(tmplm.lower,[],length(t));
else
    tmplm.lower=[];
end
if ~isempty(ub)
    tmplm.upper=reshape(tmplm.upper,[],length(t));
else
    tmplm.upper=[];
end
if ~isempty(A)
    tmplm.ineqlin=reshape(tmplm.ineqlin,[],length(t));
else
    tmplm.ineqlin=[];
end
lm=[tmplm.lower;tmplm.upper;tmplm.ineqlin];
dv=OCGRADCONT.OPTIONS.gradientmappinggamma*(xq-v);


function [Y_best,v_best,impr,lineiter]=linesearchloc(t,Y,v,dv,par,R,cfv,norm_dv)
global OCGRADCONT OCGRADSOL

% Y consists of states, costates and objective path
impr=false;
T=t(OCGRADCONT.TIMEMESH.num);
v0=v;
mu=OCGRADCONT.OPTIONS.initlinestepwidth;
oval0=objectivevalue(T,Y,par,v0);

phi0=oval0-R*cfv;
beta=OCGRADCONT.OPTIONS.gamma*norm_dv.^2;


lineiter=1;
while ~impr && lineiter<=OCGRADCONT.OPTIONS.maxlinesearchiter

    step=mu^lineiter;
    v_new=v0+step*dv;
    
    Y_new=OCGRADCONT.gradientodesolver(1,OCGRADSOL.statedynamicsext,t,Y(:,1),v_new,par);
    oval_new=objectivevalue(T,Y_new,par,v_new);

    if OCGRADCONT.cconstraint
        cstr=OCGRADSOL.constraint(t,Y_new,v_new,par);
        cstr(cstr>0)=0;
    else
        cstr=0;
    end
    cfv_new=max(0,abs(min(cstr(:))));

    phi_new=oval_new-R*cfv_new;
    if phi_new>phi0+step*beta
        impr=true;
        v_best=v_new;
        Y_best=Y_new;
    end
    lineiter=lineiter+1;
end
if ~impr
    % If there is no improvement, return original values
    v_best=v0;
    Y_best=Y;
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
    cst_y=OCGRADSOL.transversaltycondition(T,Y(OCGRADCONT.statecoordinate.concentrated,OCGRADCONT.TIMEMESH.num),par);
end

