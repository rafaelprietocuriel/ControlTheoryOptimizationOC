function out=extremalgrad4ft()

out{1}=@operatoreq;
out{2}=@frechetder;
out{3}=@gradientsolution;
out{7}=@defaultprocessor;
out{8}=@testfunc;
out{9}=@targetvaluefunc;
out{10}=@probleminit;
out{11}=@operatorpfrechet;
out{16}=@done;
out{17}=@adapt;
out{18}=@predictextremal;
out{19}=@predictextremaldiff;
out{20}=@workspaceadapt;
out{21}=@plotcontinuation;
out{22}=@formatsolution;
out{23}=@testadmissibility;
out{24}=@drearr;
out{25}=@plotsolution;
out{26}=@saveintermediate;
out{27}=@datapath;
out{30}=@printcontinuation;

function [res,extremal,graditer,dgraditer,lineiter,dlineiter,numrepdd]=operatoreq(coeff,extremal,tangent)
global OCGRADCONT OCGRADSOL

res=[];
graditer=[];
dgraditer=[];
lineiter=[];
dlineiter=[];
numrepdd=[];
switch OCGRADCONT.conttype
    case 'endstate'
        [freepar,modelpar]=drearr(coeff);
        yT=OCGRADCONT.endstate;

        yT(OCGRADCONT.fixendstatecoordinate)=yT(OCGRADCONT.fixendstatecoordinate)+freepar(OCGRADCONT.continuationindex)*OCGRADCONT.targetvector; % adapt end state values
    otherwise
        return
end

extremal.cst_y(OCGRADCONT.freeendcostatecoordinate,end)=freepar(OCGRADCONT.freeendcostateindex); % set final costate values to the free parameter values

[extremal,graditer,lineiter]=OCGRADCONT.gradientsolver(extremal,modelpar);

[extremal.tangent,dgraditer,dlineiter,numrepdd]=directionalderivative(coeff,extremal);

% condition at the final time
res=yT-extremal.y((OCGRADCONT.fixendstatecoordinate),end);
res(OCGRADCONT.HE.numdvariables,1)=0;


function [res,extremal,graditer,dgraditer,lineiter,dlineiter,numrepdd]=operatoreqred(coeff,extremal,tangent)
global OCGRADCONT OCGRADSOL

res=[];
graditer=[];
dgraditer=[];
lineiter=[];
dlineiter=[];
numrepdd=[];
switch OCGRADCONT.conttype
    case 'endstate'
        [freepar,modelpar]=drearr(coeff);
        yT=OCGRADCONT.endstate;

        yT(OCGRADCONT.fixendstatecoordinate)=yT(OCGRADCONT.fixendstatecoordinate)+freepar(OCGRADCONT.continuationindex)*OCGRADCONT.targetvector; % adapt end state values
    otherwise
        return
end

extremal.cst_y(OCGRADCONT.freeendcostatecoordinate,end)=freepar(OCGRADCONT.freeendcostateindex); % set final costate values to the free parameter values

[extremal,graditer,lineiter]=OCGRADCONT.gradientsolver(extremal,modelpar);

% condition at the final time
res=yT-extremal.y((OCGRADCONT.fixendstatecoordinate),end);

function varargout=probleminit(varargin)
coeff=varargin{1};
tangent=varargin{2};
WorkspaceInit(coeff,tangent);

% all done succesfully
varargout{1} = 0;

function J=frechetder(coeff,extremal,tangent)
global OCGRADCONT OCGRADSOL
J=[];
switch OCGRADCONT.conttype
    case 'endstate'
        [freepar,modelpar]=drearr(coeff);
        yT=OCGRADCONT.endstate;

        yT(OCGRADCONT.fixendstatecoordinate)=yT(OCGRADCONT.fixendstatecoordinate)+freepar(OCGRADCONT.continuationindex)*OCGRADCONT.targetvector; % adapt end state values
    otherwise
        return
end
%DYDp=reshape(extremal.tangent.y(:,OCGRADCONT.TIMEMESH.num),[],length(OCGRADCONT.fixendstatecoordinate));
J=OCGRADCONT.Jacobian4Equation;

%J(OCGRADCONT.fixendstatecoordinate,OCGRADCONT.fixendstatecoordinate)=DYDp;
J(OCGRADCONT.freeendcostateindex,OCGRADCONT.continuationindex)=OCGRADCONT.targetvector;
res0=operatoreqred(coeff,extremal,tangent);
% J=OCGRADCONT.Jacobian4Equation;
directionalderivativestep=OCGRADCONT.OPTIONS.directionalderivativestep;
loopidx=OCGRADCONT.freeendcostateindex;
for ii=loopidx % loop for the different solutions
    coeff1=coeff;
    coeff1(ii)=coeff1(ii)+directionalderivativestep;
    res=operatoreqred(coeff1,extremal,tangent);
    J(:,ii)=(res-res0)/directionalderivativestep;
end
%-------------------------------------------------------------------------
function out=operatorpfrechet(coeff,extremal,tangent)
out=[];

%-------------------------------------------------------------------------
function [extremal,graditer,dgraditer,lineiter,dlineiter,numrepdd,statinfo]=gradientsolution(coeff,extremal,tangent)
global OCGRADCONT OCGRADSOL
[freepar,modelpar]=drearr(coeff);
switch OCGRADCONT.conttype
    case 'time'
        if OCGRADCONT.contindex==1
        elseif OCGRADCONT.contindex==2
            extremal.t=extremal.t*freepar(OCGRADCONT.continuationindex);
            extremal.timehorizon=freepar;
        end
    case 'initialstate'
        y0=OCGRADCONT.initialstate;
        y0(OCGRADCONT.initialstatecoordinate)=y0(OCGRADCONT.initialstatecoordinate)+freepar(OCGRADCONT.continuationindex)*OCGRADCONT.targetvector;
        extremal.y=y0;
end
[extremal,graditer,lineiter,statinfo]=OCGRADCONT.gradientsolver(extremal,modelpar);
switch OCGRADCONT.conttype
    case 'time'
        extremal.t=extremal.t/extremal.t(end);
end
[extremal.tangent,dgraditer,dlineiter,numrepdd]=directionalderivative(coeff,extremal);
extremal.coeff=[extremal.v(:);extremal.y(:,1);extremal.cst_y(:,end)];
switch OCGRADCONT.conttype
    case 'time'
        extremal.t=extremal.t/extremal.t(end);
end

%-------------------------------------------------------------------------
function [precoeff,extremal]=predictextremal(coeff,extremal,stepwidth,tangent)
global OCGRADCONT OCGRADSOL
tottangent=[sign(tangent(end))*extremal.tangent.coeff(:,end);tangent];
tottangent=tottangent/norm(tottangent);
tot_coeff=[extremal.coeff;coeff];
tot_coeff=tot_coeff+stepwidth*tottangent;
precoeff=tot_coeff((end-OCGRADCONT.HE.numdvariables+1):end);
extremal.y=tot_coeff(OCGRADCONT.initiallocstatecoordinate);
extremal.cst_y=tot_coeff(OCGRADCONT.endloccostatecoordinate);
if OCGRADCONT.control_num.concentrated==1
    extremal.v=tot_coeff(OCGRADCONT.loccontrolcoordinate).';
else
    extremal.v=tot_coeff(OCGRADCONT.loccontrolcoordinate);
end
extremal.tangent=[];

%-------------------------------------------------------------------------
function [coeffpre,extremal]=predictextremaldiff(coeff1,coeff2,extremal1,extremal2,fac)
global OCGRADCONT OCGRADSOL
coeffpre=coeff1+fac*(coeff2-coeff1);
freepar=drearr(coeffpre);
extremal=extremal1;
switch OCGRADCONT.conttype
    case 'time'
        if OCGRADCONT.contindex==1
        elseif OCGRADCONT.contindex==2
            extremal.t=extremal.t*freepar;
            extremal.timehorizon=freepar;
        end
end

extremal.y=extremal1.y+fac*(extremal2.y-extremal1.y);
extremal.v=extremal1.v+fac*(extremal2.v-extremal1.v);
extremal.cst_y=extremal1.cst_y+fac*(extremal2.cst_y-extremal1.cst_y);
%[extremal,graditer,totallineiter]=OCGRADCONT.gradientsolver(extremal,modelpar);
switch OCGRADCONT.conttype
    case 'time'
        extremal.t=extremal.t/extremal.t(end);
end
%----------------------------------------------------------------
function h=plotcontinuation(coeff,extremal,tangent,contdata,makemovie)
global OCGRADSOL OCGRADCONT
[freepar,modelpar]=drearr(coeff);
switch OCGRADCONT.conttype
    case 'time'
        if OCGRADCONT.contindex==1
        elseif OCGRADCONT.contindex==2
            extremal.t=extremal.t*freepar;
        end
end
h=OCGRADSOL.plotcontinuation(freepar,extremal,tangent,modelpar,contdata,OCGRADCONT.conttype,makemovie);

%----------------------------------------------------------------
function h=plotsolution(coeff,extremal,tangent)
global OCGRADSOL OCGRADCONT
[freepar,modelpar]=drearr(coeff);
h=OCGRADSOL.plotsolution(coeff,extremal,tangent,modelpar);


%----------------------------------------------------------------
function idx=printcontinuation(coeff,extremal,tangent)
global OCGRADCONT OCGRADSOL
idx=[];
if isempty(coeff)
    return
end
if isempty(OCGRADCONT.targetvalue)
    fprintf(1,' Continuation parameter: %g\n',coeff(OCGRADCONT.continuationindex));
else
    out=targetvaluefunc(1,coeff,extremal,tangent);
    fprintf(1,' Distance from Target value: %g\n',out);
end

%---------------------------------------------------------------
function varargout=defaultprocessor(varargin)
global OCGRADCONT OCGRADSOL

varargout{2}=nan;
% all done succesfully
varargout{1}=0;
%-------------------------------------------------------------
function [out, failed]=testfunc(id,coeff,extremal,tangent)
global OCGRADCONT OCGRADSOL

out(1)=0;
failed=[];
%----------------------------------------------------------------
function [out, failed]=targetvaluefunc(id,coeff,extremal,tangent)
global OCGRADCONT OCGRADSOL

failed=[];
for ii=id
    switch ii
        case 1
            if isempty(OCGRADCONT.targetvalue)
                out=[];
            else
                switch OCGRADCONT.conttype
                    case 'time'
                        freepar=drearr(coeff);
                        out=OCGRADCONT.targetvalue-freepar;
                    case 'initialstate'
                        freepar=drearr(coeff);

                        out=1-freepar;
                    case 'endstate'
                        freepar=drearr(coeff);

                        out=1-freepar(end);
                    case 'parameter'
                        freepar=drearr(coeff);

                        out=OCGRADCONT.targetvalue-freepar;
                end
            end

        otherwise
            out=[];
            failed=1;
    end
end

%----------------------------------------------------------------
function [b,info,label]=testadmissibility(coeff,extremal,tangent)
global OCGRADCONT OCGRADSOL 
info=[];
label=[];

[freepar,modelpar]=drearr(coeff);
switch OCGRADCONT.conttype
    case 'time'
        if OCGRADCONT.contindex==1
        elseif OCGRADCONT.contindex==2
            extremal.t=extremal.t*freepar;
            extremal.timehorizon=freepar;
        end
end
val=OCGRADSOL.admissible(extremal.t,extremal.y,extremal.v,modelpar);
b=any(val(:)<-OCGRADCONT.OPTIONS.zerodeviationtolerance);

%----------------------------------------------------------------
function out=formatsolution(coeff,extremal,tangent)
global OCGRADCONT OCGRADSOL 
[freepar,modelpar]=drearr(coeff);
switch OCGRADCONT.conttype
    case 'time'
        if OCGRADCONT.contindex==1
        elseif OCGRADCONT.contindex==2
            extremal.timehorizon=freepar;
            extremal.t=extremal.t*freepar;
        end
end

out.extremal=extremal;
out.parameters=freepar;
out.extremal.modelparameter=modelpar;
out.extremal.modelname=OCGRADCONT.modelname;
out.extremal.type=OCGRADCONT.octype;
out.extremal.solver.gradient=OCGRADCONT.gradientsolver;
out.extremal.solver.ode=OCGRADCONT.gradientodesolver;
out.extremal.solver.option.R0=OCGRADCONT.OPTIONS.R0;
out.extremal.solver.option.gradientmappingmethod=OCGRADCONT.OPTIONS.gradientmappingmethod;
out.extremal.solver.option.gradtol=OCGRADCONT.OPTIONS.gradtol;
out.extremal.solver.option.gamma=OCGRADCONT.OPTIONS.gamma;
out.extremal.solver.option.mu=OCGRADCONT.OPTIONS.mu;
out.extremal.contclass=OCGRADCONT.problem_func;
out.extremal.conttype=OCGRADCONT.conttype;
switch OCGRADCONT.conttype
    case 'time'
    case 'initialstate'
        out.extremal.initialstate=OCGRADCONT.initialstate;
        out.extremal.targetvector=OCGRADCONT.targetvector;
        out.extremal.initialstatecoordinate=OCGRADCONT.initialstatecoordinate;
        out.extremal.continuationindex=OCGRADCONT.continuationindex;
    case 'parameter'
        out.extremal.freeparameterindex=OCGRADCONT.contindex;
end
%---------------------------------------------------------
function done
WorkspaceDone;
% ---------------------------------------------------------

function WorkspaceInit(coeff,extremal,tangent)
global OCGRADCONT OCGRADSOL

% ------------------------------------------------------

function WorkspaceDone

% ------------------------------------------------------

function [freepar,modelpar]=drearr(coeff)
global OCGRADCONT OCGRADSOL

modelpar=OCGRADCONT.modelparameter;
freepar=coeff;
switch  OCGRADCONT.conttype
    case 'parameter'
        modelpar(OCGRADCONT.contindex)=freepar;%OCGRADCONT.startingparametervalue+freepar*OCGRADCONT.targetvector;
end

%-----------------------------------------------------------------
function failed=saveintermediate(sout,gradout,contnum)
global OCGRADCONT OCGRADSOL
failed=0;
MODELINFO.OCGRADCONT=OCGRADCONT;
MODELINFO.OCGRADSOL=OCGRADSOL;
try
    if contnum==1
        save([OCGRADSOL.basicglobalvarfilename '4extremalgrad4ft'],'MODELINFO')
    end
    save([OCGRADSOL.basicresultfilename '4extremalgrad4ft'],'sout','gradout')
catch
    failed=1;
end

%-----------------------------------------------------------------
function pathname=datapath()
global OCGRADSOL

pathname=OCGRADSOL.datapath();

% helper function for directional derivative
function [dX,graditer,lineiter,ctr]=directionalderivative(coeff,extremal)
global OCGRADCONT

directionalderivativestep=OCGRADCONT.OPTIONS.directionalderivativestep;
ctr=1;
%linestepwidth0=OCGRADCONT.linestepwidth;
%OCGRADCONT.linestepwidth=linestepwidth0*1e-3;
y=extremal.y;
v=extremal.v;
cst_y=extremal.cst_y;
max_graditer=-inf;
max_lineiter=-inf;
dX.y=[];
dX.v=[];
dX.cst_y=[];
dX.coeff=[];
for ii=1:length(coeff) % loop for the different solutions
    coeff1=coeff;
    coeff1(ii)=coeff1(ii)+directionalderivativestep;

    [freepar,modelpar]=drearr(coeff1);
    switch OCGRADCONT.conttype
        case 'time'
            if OCGRADCONT.contindex==1
            elseif OCGRADCONT.contindex==2
                extremal.timehorizon=freepar;
                extremal.t=extremal.t*freepar;
            end
        case 'initialstate'
            y0=OCGRADCONT.initialstate;
            y0(OCGRADCONT.initialstatecoordinate)=y0(OCGRADCONT.initialstatecoordinate)+freepar(OCGRADCONT.continuationindex)*OCGRADCONT.targetvector;
            extremal.y=y0;
    case 'endstate'
        extremal.cst_y(OCGRADCONT.freeendcostatecoordinate,end)=freepar(OCGRADCONT.freeendcostateindex); % set final costate values to the free parameter values

    end
    [extremal1,graditer,lineiter]=OCGRADCONT.gradientsolver(extremal,modelpar);
    switch OCGRADCONT.conttype
        case 'time'
            extremal.t=linspace(0,1,length(extremal.t));
    end
    dx.y=(extremal1.y-y)/directionalderivativestep;
    dx.v=(extremal1.v-v)/directionalderivativestep;
    dx.cst_y=(extremal1.cst_y-cst_y)/directionalderivativestep;
    dx.coeff=[dx.v(:);dx.y(:,1);dx.cst_y(:,OCGRADCONT.TIMEMESH.num)];
    max_graditer=max(max_graditer,graditer);
    max_lineiter=max(max_lineiter,lineiter);
    dX.y=[dX.y;dx.y];
    dX.cst_y=[dX.cst_y;dx.cst_y];
    dX.v=[dX.v;dx.v];
    dX.coeff=[dX.coeff,dx.coeff];
end
%OCGRADCONT.linestepwidth=linestepwidth0;

