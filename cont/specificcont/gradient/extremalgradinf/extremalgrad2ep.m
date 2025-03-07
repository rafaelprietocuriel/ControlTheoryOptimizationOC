function out=extremalgrad2ep()

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

function [res,extremal,graditer,dgraditer,lineiter,dlineiter,numrepdd,statinfo]=operatoreq(coeff,extremal,tangent)
global OCGRADCONT OCGRADSOL
[freepar,modelpar]=drearr(coeff);
y0=OCGRADCONT.initialstatevalue;
switch OCGRADCONT.conttype
    case 'state'
        y0=y0+freepar(OCGRADCONT.continuationcoordinate)*OCGRADCONT.contvector; % adapt initial state values
end
if OCGRADCONT.normalizetime
    extremal.t=extremal.t*freepar(OCGRADCONT.freetimehorizoncoordinate);
end
extremal.y=y0; % set initial state values
extremal.cst_y(:,end)=freepar(OCGRADCONT.freeparametercoordinate); % set final costate values to the free parameter values
[extremal,graditer,lineiter,statinfo]=OCGRADCONT.gradientsolver(extremal,modelpar);
if OCGRADCONT.normalizetime
    extremal.t=extremal.t/extremal.t(OCGRADCONT.TIMEMESH.num);
end

[extremal.tangent,dgraditer,dlineiter,numrepdd]=directionalderivative(coeff,extremal);

% asymptotic condition at the final time
res=OCGRADSOL.asymptoticmatrix.'*(OCGRADSOL.saddlepoint-[extremal.y(:,OCGRADCONT.TIMEMESH.num);freepar(OCGRADCONT.freeparametercoordinate)]);
if OCGRADCONT.dist2ep
    % keep the solution in a constant distance from the equilibrium
    res(OCGRADCONT.freetimehorizoncoordinate,1)=OCGRADCONT.dist2epvalue-sqrt(sum((OCGRADSOL.saddlepoint-[extremal.y(:,OCGRADCONT.TIMEMESH.num);freepar(OCGRADCONT.freeparametercoordinate)]).^2));
end
res(OCGRADCONT.HE.numdvariables,1)=0;

%-------------------------------------------------------------------------
function [extremal,graditer,dgraditer,lineiter,dlineiter,numrepdd]=gradientsolution(coeff,extremal,tangent)
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
[extremal,graditer,lineiter]=OCGRADCONT.gradientsolver(extremal,modelpar);
switch OCGRADCONT.conttype
    case 'time'
        extremal.t=linspace(0,1,length(extremal.t));
end
[extremal.tangent,dgraditer,dlineiter,numrepdd]=directionalderivative(coeff,extremal);
extremal.coeff=[extremal.v(:);extremal.y(:,1);extremal.cst_y(:,end)];
switch OCGRADCONT.conttype
    case 'time'
        extremal.t=linspace(0,1,length(extremal.t));
end

function res=operatoreqred(coeff,extremal,tangent)
global OCGRADCONT OCGRADSOL
[freepar,modelpar]=drearr(coeff);
y0=OCGRADCONT.initialstatevalue;
switch OCGRADCONT.conttype
    case 'state'
        y0=y0+freepar(OCGRADCONT.continuationcoordinate)*OCGRADCONT.contvector;
end
extremal.y=y0;
extremal.cst_y(:,end)=freepar(OCGRADCONT.freeparametercoordinate);

if OCGRADCONT.normalizetime
    extremal.t=extremal.t*freepar(OCGRADCONT.freetimehorizoncoordinate);
end
extremal=OCGRADCONT.gradientsolver(extremal,modelpar);
if OCGRADCONT.normalizetime
    extremal.t=extremal.t/extremal.t(OCGRADCONT.TIMEMESH.num);
end

res=OCGRADSOL.asymptoticmatrix.'*(OCGRADSOL.saddlepoint-[extremal.y(:,end);freepar(OCGRADCONT.freeparametercoordinate)]);
if OCGRADCONT.dist2ep
    % keep the solution in a constant distance from the equilibrium
    res(OCGRADCONT.freetimehorizoncoordinate,1)=OCGRADCONT.dist2epvalue-sqrt(sum((OCGRADSOL.saddlepoint-[extremal.y(:,OCGRADCONT.TIMEMESH.num);freepar(OCGRADCONT.freeparametercoordinate)]).^2));
end


function varargout=probleminit(varargin)
coeff=varargin{1};
tangent=varargin{2};
WorkspaceInit(coeff,tangent);

% all done succesfully
varargout{1} = 0;

function J=frechetder(coeff,extremal,tangent)
global OCGRADCONT OCGRADSOL
[freepar,modelpar]=drearr(coeff);

J=OCGRADCONT.Jacobian4Equation;

J(OCGRADCONT.freeparametercoordinate,OCGRADCONT.freeparametercoordinate)=-OCGRADSOL.asymptoticmatrix(OCGRADCONT.costatecoordinate.concentrated,OCGRADCONT.statecoordinate.concentrated).';
%J0=J;
DYDp=[reshape(extremal.tangent.y(:,OCGRADCONT.TIMEMESH.num),[],OCGRADCONT.continuationcoordinate); ...
             reshape(extremal.tangent.cst_y(:,OCGRADCONT.TIMEMESH.num),[],OCGRADCONT.continuationcoordinate)];
J(OCGRADCONT.freeparametercoordinate,OCGRADCONT.continuationcoordinate)=-OCGRADSOL.asymptoticmatrix.'*DYDp(:,OCGRADCONT.continuationcoordinate);         
if OCGRADCONT.dist2ep
    DXT=OCGRADSOL.saddlepoint-[extremal.y(:,OCGRADCONT.TIMEMESH.num);freepar(OCGRADCONT.freeparametercoordinate)];
    J(OCGRADCONT.freeparametercoordinate,OCGRADCONT.freetimehorizoncoordinate)=-OCGRADSOL.asymptoticmatrix.'*DYDp(:,OCGRADCONT.freetimehorizoncoordinate);         
    J(OCGRADCONT.freetimehorizoncoordinate,:)=(DXT.'*DYDp/(sqrt(sum(DXT.^2))));
end
% res0=operatoreqred(coeff,extremal,tangent);
% %J0=J;
% if OCGRADCONT.dist2ep
%     %loopidx=(OCGRADCONT.continuationcoordinate-1):OCGRADCONT.continuationcoordinate;
%     loopidx=1:OCGRADCONT.continuationcoordinate;
% else
%     loopidx=1:OCGRADCONT.continuationcoordinate;
% end
% directionalderivativestep=OCGRADCONT.OPTIONS.directionalderivativestep;
% for ii=loopidx % loop for the different solutions
%     coeff1=coeff;
%     coeff1(ii)=coeff1(ii)+directionalderivativestep;
%     res=operatoreqred(coeff1,extremal,tangent);
%     J(:,ii)=(res-res0)/directionalderivativestep;
% end
% J=J0;
 %dJ=J-J0;
% if max(abs(dJ(:)))>1e-7
%     J-J0
% end

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
tangent=tottangent(end-length(coeff)+1:end);
%-------------------------------------------------------------------------
function [coeffpre,extremal]=predictextremaldiff(coeff1,coeff2,extremal1,extremal2,fac)
global OCGRADCONT OCGRADSOL
coeffpre=coeff1+fac*(coeff2-coeff1);

extremal=extremal1;

extremal.y=extremal1.y+fac*(extremal2.y-extremal1.y);
extremal.v=extremal1.v+fac*(extremal2.v-extremal1.v);
extremal.cst_y=extremal1.cst_y+fac*(extremal2.cst_y-extremal1.cst_y);
%extremal=OCGRADCONT.gradientsolution(coeffpre,extremal,[]);
%----------------------------------------------------------------
function h=plotcontinuation(coeff,extremal,tangent,contdata,makemovie)
global OCGRADSOL OCGRADCONT
[freepar,modelpar]=drearr(coeff);
if OCGRADCONT.normalizetime
    extremal.t=extremal.t*freepar(OCGRADCONT.freetimehorizoncoordinate);
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
    fprintf(1,' Continuation parameter: %g\n',coeff(OCGRADCONT.continuationcoordinate));
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
                freepar=drearr(coeff);
                switch OCGRADCONT.conttype
                    case 'state'
                        out=1-freepar(OCGRADCONT.continuationcoordinate);
                    case 'time'
                        switch OCGRADCONT.targettype
                            case 'contpar'
                                out=OCGRADCONT.targetvalue-freepar(OCGRADCONT.continuationcoordinate);
                            case 'norm2eq'
                                out=OCGRADCONT.targetvalue-sqrt(sum((OCGRADSOL.saddlepoint-[extremal.y(:,OCGRADCONT.TIMEMESH.num);extremal.cst_y(:,OCGRADCONT.TIMEMESH.num)]).^2));
                        end
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

bdist=0;
b=0;

[freepar,modelpar]=drearr(coeff);
switch OCGRADCONT.conttype
    case 'time'
        if OCGRADCONT.contindex==1
        elseif OCGRADCONT.contindex==2
            extremal.t=extremal.t*freepar(OCGRADCONT.continuationcoordinate);
            extremal.timehorizon=freepar(OCGRADCONT.continuationcoordinate);
        end
end
val=OCGRADSOL.admissible(extremal.t,extremal.y,extremal.v,modelpar);
if OCGRADCONT.testdistance
    bdist=(OCGRADCONT.OPTIONS.maxdistance-sqrt(sum((OCGRADSOL.saddlepoint-[extremal.y(:,OCGRADCONT.TIMEMESH.num);freepar(OCGRADCONT.freeparametercoordinate)]).^2)))<0;
end
if bdist
    label='dist';
    b=bdist;
elseif any(val(:)<-OCGRADCONT.OPTIONS.zerodeviationtolerance)
    b=1;
    label='cstr';
end

%----------------------------------------------------------------
function out=formatsolution(coeff,extremal,tangent)
global OCGRADCONT OCGRADSOL
[freepar,modelpar]=drearr(coeff);

if OCGRADCONT.normalizetime
    %extremal.t=extremal.t/extremal.t(end);
    extremal.t=extremal.t*freepar(OCGRADCONT.freetimehorizoncoordinate);
    extremal.timehorizon=freepar(OCGRADCONT.freetimehorizoncoordinate);
end

out.extremal=extremal;
out.extremal.freeparparameter=freepar;
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
out.extremal.equilibrium=OCGRADSOL.saddlepoint;
%out.extremal.arcarg4equilibrium=OCGRADSOL.arcarg4equilibrium;
out.extremal.freeparametercoordinate=OCGRADCONT.freeparametercoordinate;
out.extremal.continuationcoordinate=OCGRADCONT.continuationcoordinate;
out.extremal.dist2ep=OCGRADCONT.dist2ep;
if OCGRADCONT.dist2ep
    out.extremal.freetimehorizoncoordinate=OCGRADCONT.freetimehorizoncoordinate;
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

%-----------------------------------------------------------------
function failed=saveintermediate(sout,gradout,contnum)
global OCGRADCONT OCGRADSOL
failed=0;
MODELINFO.OCGRADCONT=OCGRADCONT;
MODELINFO.OCGRADSOL=OCGRADSOL;
try
    if contnum==1
        save([OCGRADSOL.basicglobalvarfilename '4extremalgrad2ep'],'MODELINFO')
    end
    save([OCGRADSOL.basicresultfilename '4extremalgrad2ep'],'sout','gradout')
catch
    failed=1;
end

%-----------------------------------------------------------------
function pathname=datapath()
global OCGRADSOL

pathname=OCGRADSOL.datapath();

% helper function for directional derivative
function [dX,max_graditer,max_lineiter,ctr_max]=directionalderivative(coeff,extremal)
global OCGRADCONT

extremal0=extremal;
ctrder=0;
y=extremal.y;
v=extremal.v;
cst_y=extremal.cst_y;
o=extremal.o;
max_graditer=-inf;
max_lineiter=-inf;
ctr_max=-inf;
dX.y=[];
dX.v=[];
dX.cst_y=[];
dX.o=[];
dX.coeff=[];
% for ii=1:length(coeff) % loop for the different solutions
%     extremal=extremal0;
%     ctrder=ctrder+1;
%     directionalderivativestep=OCGRADCONT.OPTIONS.directionalderivativestep;
%     ctr=0;
%     while 1
%         ctr=ctr+1;
%         coeff1=coeff;
%         coeff1(ii)=coeff1(ii)+directionalderivativestep;
%         [freepar,modelpar]=drearr(coeff1);
%         y0=OCGRADCONT.initialstatevalue;
%         switch OCGRADCONT.conttype
%             case 'state'
%                 y0=y0+freepar(OCGRADCONT.continuationcoordinate)*OCGRADCONT.contvector;
%         end
%         extremal.y=y0;
%         extremal.cst_y(:,end)=freepar(OCGRADCONT.freeparametercoordinate);
%         if OCGRADCONT.normalizetime
%             extremal.t=extremal.t*freepar(OCGRADCONT.freetimehorizoncoordinate);
%         end
%         [extremal1,graditer,lineiter]=OCGRADCONT.gradientsolver(extremal,modelpar);
%         dx.y=(extremal1.y-y)/directionalderivativestep;
%         dx.v=(extremal1.v-v)/directionalderivativestep;
%         dx.cst_y=(extremal1.cst_y-cst_y)/directionalderivativestep;
%         dx.o=(extremal1.o-o)/directionalderivativestep;
%         dx.coeff=[dx.v(:);dx.y(:,1);dx.cst_y(:,OCGRADCONT.TIMEMESH.num)];
%         if norm(dx.coeff)~=0
%             max_graditer=max(max_graditer,graditer);
%             max_lineiter=max(max_lineiter,lineiter);
%             ctr_max=max(ctr_max,ctr);
%             dX.y=[dX.y;dx.y];
%             dX.cst_y=[dX.cst_y;dx.cst_y];
%             dX.v=[dX.v;dx.v];
%             dX.o=[dX.o;dx.o];
%             dX.coeff=[dX.coeff,dx.coeff];
%             break
%         end
%         if ctr>2
%             break
%         end
%         directionalderivativestep=directionalderivativestep*5;
%         extremal=extremal0;
%     end
% end
for ii=1:length(coeff) % loop for the different solutions
    extremal=extremal0;
    ctrder=ctrder+1;
    directionalderivativestep=OCGRADCONT.OPTIONS.directionalderivativestep;
    coeff1=coeff;
    coeff1(ii)=coeff1(ii)+directionalderivativestep;
    [freepar,modelpar]=drearr(coeff1);
    y0=OCGRADCONT.initialstatevalue;
    switch OCGRADCONT.conttype
        case 'state'
            y0=y0+freepar(OCGRADCONT.continuationcoordinate)*OCGRADCONT.contvector;
    end
    extremal.y=y0;
    extremal.cst_y(:,end)=freepar(OCGRADCONT.freeparametercoordinate);
    if OCGRADCONT.normalizetime
        extremal.t=extremal.t*freepar(OCGRADCONT.freetimehorizoncoordinate);
    end
    [extremal1,graditer,lineiter]=OCGRADCONT.gradientsolver(extremal,modelpar);
    dx.y=(extremal1.y-y)/directionalderivativestep;
    dx.v=(extremal1.v-v)/directionalderivativestep;
    dx.cst_y=(extremal1.cst_y-cst_y)/directionalderivativestep;
    dx.o=(extremal1.o-o)/directionalderivativestep;
    dx.coeff=[dx.v(:);dx.y(:,1);dx.cst_y(:,OCGRADCONT.TIMEMESH.num)];

    max_graditer=max(max_graditer,graditer);
    max_lineiter=max(max_lineiter,lineiter);
    dX.y=[dX.y;dx.y];
    dX.cst_y=[dX.cst_y;dx.cst_y];
    dX.v=[dX.v;dx.v];
    dX.o=[dX.o;dx.o];
    dX.coeff=[dX.coeff,dx.coeff];
end