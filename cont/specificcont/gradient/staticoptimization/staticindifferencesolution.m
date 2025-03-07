function out=staticindifferencesolution()

out{1}=@operatoreq;
out{2}=@frechetder;
out{3}=@gradientsolution;
out{7}=@defaultprocessor;
out{8}=@testfunc;
out{9}=@targetvaluefunc;
out{10}=@probleminit;
out{11}=@operatorpfrechet;
out{13}=@singmat;
out{14}=@process;
out{15}=@locate;
out{16}=@done;
out{17}=@adapt;
out{18}=@predictextremal;
out{19}=@predictextremaldiff;
out{20}=@workspaceadapt;
out{21}=@plotcontinuation;
out{22}=@formatsolution;
out{24}=@drearr;
out{25}=@plotsolution;
out{26}=@saveintermediate;
out{27}=@datapath;
out{30}=@printcontinuation;

function [res,extremal,totgraditer]=operatoreq(coeff,extremal,tangent)
global OCGRADCONT OCGRADSOL
res=zeros(OCGRADSOL.degree,1);

[freepar,modelpar]=drearr(coeff);

totgraditer=0;
for ii=1:OCGRADSOL.degree
    [extremal(ii),graditer]=OCGRADCONT.gradientsolver(extremal(ii),modelpar);
    totgraditer=totgraditer+graditer;
end
for ii=1:OCGRADSOL.degree % last parameter is the continuaton parameter
    extremal(ii).tangent=directionalderivative(coeff,extremal(ii));
    if ii<OCGRADSOL.degree
        res(ii)=extremal(ii+1).objectivevalue-extremal(ii).objectivevalue;
    end
end


function res=operatoreqred(coeff,extremal,tangent)
global OCGRADCONT OCGRADSOL
res=zeros(OCGRADSOL.degree-1,1);

[freepar,modelpar]=drearr(coeff);

for ii=1:OCGRADSOL.degree
    extremal(ii)=OCGRADCONT.gradientsolver(extremal(ii),modelpar);
    if ii>1
        res(ii-1)=extremal(ii).objectivevalue-extremal(ii-1).objectivevalue;
    end
end

function J=frechetder(coeff,extremal,tangent,varargin)
global OCGRADCONT OCGRADSOL
% this step doubles 'operatoreq', better solution to be implemented
% [freepar,modelpar]=drearr(coeff);
% for ii=1:OCGRADSOL.degree
%     gradobj(:,ii)=OCGRADSOL.gradientobjective(extremal(ii).y,modelpar);
%     pargradobj(:,ii)=OCGRADSOL.gradientparobjective(extremal(ii).y,modelpar);
%     dx(:,:,ii)=directionalderivative(coeff,extremal(ii));
% end
% 
% pargradobj=pargradobj([OCGRADCONT.freeparameterindex OCGRADCONT.contindex],:);
% 
% J=[(gradobj(:,2).'*dx(:,2,1)+pargradobj(2,1))-(gradobj(:,1).'*dx(:,1,1)+pargradobj(1,1)), ...
%     (gradobj(:,2).'*dx(:,2,2)+pargradobj(2,2))-(gradobj(:,1).'*dx(:,1,2)+pargradobj(2,1))];

J=zeros(OCGRADSOL.degree-1,OCGRADSOL.degree);
res0=operatoreqred(coeff,extremal,tangent);
for ii=1:length(coeff) % loop for the different solutions
    directionalderivativestep=OCGRADCONT.OPTIONS.directionalderivativestep;
    coeff1=coeff;
    coeff1(ii)=coeff1(ii)+directionalderivativestep;
    res=operatoreqred(coeff1,extremal,tangent);
    J(1,ii)=(res-res0)/directionalderivativestep;
end

function varargout=probleminit(varargin)
coeff=varargin{1};
tangent=varargin{2};
WorkspaceInit(coeff,tangent);

% all done succesfully
varargout{1} = 0;


%-------------------------------------------------------------------------
function [precoeff,extremal]=predictextremal(coeff,extremal,stepwidth,tangent)
global OCGRADCONT OCGRADSOL
tot_extremal_tangent=[];
tot_extremal_coeff=[];
for ii=1:OCGRADSOL.degree
    tot_extremal_tangent=[tot_extremal_tangent;extremal(ii).tangent(:,OCGRADSOL.degree)];
    tot_extremal_coeff=[tot_extremal_coeff;extremal(ii).y(:)];
end
tot_extremal_tangent=tot_extremal_tangent*tangent(end);
tot_extremal_tangent=[tot_extremal_tangent;tangent];
tot_extremal_tangent=tot_extremal_tangent/norm(tot_extremal_tangent);
tot_extremal_coeff=[tot_extremal_coeff;coeff];
tot_precoeff=tot_extremal_coeff+stepwidth*tot_extremal_tangent;
precoeff=tot_precoeff((end-OCGRADSOL.degree+1):end);
for ii=1:OCGRADSOL.degree
    extremal(ii).y=tot_precoeff(OCGRADCONT.statecoordinate{ii});
end

%-------------------------------------------------------------------------
function [coeffpre,extremal]=predictextremaldiff(coeff1,coeff2,extremal1,extremal2,fac)
global OCGRADCONT OCGRADSOL
coeffpre=coeff1+fac*(coeff2-coeff1);
extremal=extremal1;

for ii=1:OCGRADSOL.degree
    extremal(ii).y=extremal1(ii).y+fac*(extremal2(ii).y-extremal1(ii).y);
end
%----------------------------------------------------------------
function h=plotcontinuation(coeff,extremal,tangent,contdata,makemovie)
global OCGRADSOL OCGRADCONT
[freepar,modelpar]=drearr(coeff);
h=OCGRADSOL.plotcontinuation(coeff,extremal,tangent,modelpar,contdata,makemovie);


%----------------------------------------------------------------
function idx=printcontinuation(coeff,extremal,tangent)
global OCGRADCONT OCGRADSOL
idx=[];
if isempty(coeff)
    return
end
if isempty(OCGRADCONT.targetvalue)
    fprintf(1,' Continuation parameter: %g\n',coeff(OCGRADSOL.continuationcoordinate));
else
    [freepar,modelpar]=drearr(coeff);
    fprintf(1,' Distance to target value: %g\n',OCGRADCONT.targetvalue-modelpar(OCGRADCONT.targetindex));
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
                [freepar,modelpar]=drearr(coeff);
                out=OCGRADCONT.targetvalue-modelpar(OCGRADCONT.targetindex);
            end

        otherwise
            out=[];
            failed=1;
    end
end

%----------------------------------------------------------------
function out=formatsolution(coeff,extremal,tangent)
global OCGRADCONT OCGRADSOL
[freepar,modelpar]=drearr(coeff);

for ii=1:OCGRADSOL.degree
    extremal(ii).modelname=OCGRADCONT.modelname;
    extremal(ii).modelparameter=modelpar;
    extremal(ii).solver.gradient=OCGRADCONT.gradientsolver;
    extremal(ii).contclass=OCGRADCONT.problem_func;
    extremal(ii).freeparameterindex=OCGRADCONT.freeparameterindex;
    extremal(ii).contindex=OCGRADCONT.contindex;
    extremal(ii).freeparametercoordinate=OCGRADCONT.freeparametercoordinate;
    extremal(ii).continuationcoordinate=OCGRADCONT.continuationcoordinate;
    extremal(ii).freepar=freepar;
end
out.extremal=extremal;
out.parameters=freepar;
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
freepar=coeff([OCGRADCONT.freeparametercoordinate OCGRADCONT.continuationcoordinate]);

modelpar([OCGRADCONT.freeparameterindex OCGRADCONT.contindex])=freepar;

%-----------------------------------------------------------------
function failed=saveintermediate(sout,gradout,contnum)
global OCGRADCONT OCGRADSOL
failed=0;
MODELINFO.OCGRADCONT=OCGRADCONT;
MODELINFO.OCGRADSOL=OCGRADSOL;
try
    if contnum==1
        save([OCGRADSOL.basicglobalvarfilename '4staticindifferencesolution'],'MODELINFO')
    end
    save([OCGRADSOL.basicresultfilename '4staticindifferencesolution'],'sout','gradout')
catch
    failed=1;
end

%-----------------------------------------------------------------
function pathname=datapath()
global OCGRADSOL

pathname=OCGRADSOL.datapath();

% helper function for directional derivative
function dx=directionalderivative(coeff,extremal)
global OCGRADSOL OCGRADCONT

dx=zeros(OCGRADSOL.statenum,OCGRADSOL.degree);
extremal0=extremal;

for ii=1:OCGRADSOL.degree% loop for the free parameter
    coeff1=coeff;
    coeff1(ii)=coeff1(ii)+OCGRADCONT.OPTIONS.directionalderivativestep;
    [freepar,modelpar]=drearr(coeff1);
    extremal=OCGRADCONT.gradientsolver(extremal0,modelpar);
    dx(:,ii)=(extremal.y-extremal0.y)/OCGRADCONT.OPTIONS.directionalderivativestep;
end

