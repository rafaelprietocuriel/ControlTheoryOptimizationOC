function out=staticoptimization()

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

function [res,extremal,graditer]=operatoreq(coeff,extremal,tangent)
res=[];
graditer=[];

function varargout=probleminit(varargin)
coeff=varargin{1};
tangent=varargin{2};
WorkspaceInit(coeff,tangent);

% all done succesfully
varargout{1} = 0;

function J=frechetder(coeff,extremal,tangent)
J=[];

%-------------------------------------------------------------------------
function [extremal,graditer]=gradientsolution(coeff,extremal,tangent)
global OCGRADCONT OCGRADSOL
[freepar,modelpar]=drearr(coeff);
extremal.graddata.R=OCGRADSOL.actualR;
[extremal,graditer]=OCGRADCONT.gradientsolver(extremal,modelpar);
extremal.tangent=directionalderivative(coeff,extremal);

%-------------------------------------------------------------------------
function [coeff,extremal]=predictextremal(coeff,extremal,stepwidth,tangent)
global OCGRADCONT OCGRADSOL
tangent=[extremal.tangent;tangent];
tangent=tangent/norm(tangent);
totcoeff=[extremal.y;coeff];
totcoeff=totcoeff+stepwidth*tangent;
coeff=totcoeff(end);
extremal.y=totcoeff(1:OCGRADSOL.statenum);
extremal.tangent=[];


%-------------------------------------------------------------------------
function [coeffpre,extremal]=predictextremaldiff(coeff1,coeff2,extremal1,extremal2,fac)
global OCGRADCONT OCGRADSOL
coeffpre=coeff1+fac*(coeff2-coeff1);
extremal=extremal1;

extremal.y=extremal1.y+fac*(extremal2.y-extremal1.y);

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
    fprintf(1,' Continuation parameter: %g\n',coeff(OCGRADSOL.continuationindex));
else
    freepar=drearr(coeff);
    fprintf(1,' Distance to targetvalue: %g\n',OCGRADCONT.targetvalue-freepar);
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


                out=OCGRADCONT.targetvalue-freepar;
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

out.extremal=extremal;
out.freedata=freepar;
out.extremal.contclass=OCGRADCONT.problem_func;
out.extremal.freeparameterindex=OCGRADSOL.contparidx;
out.extremal.modelparameter=modelpar;
out.extremal.modelname=OCGRADSOL.modelname;
out.extremal.solver.gradientsolver=OCGRADCONT.gradientsolver;
out.extremal.solver.newtonsolver=OCGRADCONT.newtonsolver;
%---------------------------------------------------------
function done
WorkspaceDone;
% ---------------------------------------------------------

function WorkspaceInit(coeff,extremal,tangent)
global OCGRADCONT OCGRADSOL

OCGRADSOL.actualR=OCGRADCONT.OPTIONS.R0;
% ------------------------------------------------------

function WorkspaceDone

% ------------------------------------------------------

function [freepar,modelpar]=drearr(coeff)
global OCGRADCONT OCGRADSOL

modelpar=OCGRADSOL.modelparameter;
freepar=coeff(OCGRADSOL.continuationindex);

modelpar(OCGRADSOL.contparidx)=freepar;

%-----------------------------------------------------------------
function failed=saveintermediate(sout,gradout,contnum)
global OCGRADCONT OCGRADSOL
failed=0;
MODELINFO.OCGRADCONT=OCGRADCONT;
MODELINFO.OCGRADSOL=OCGRADSOL;
try
    if contnum==1
        save([OCGRADSOL.basicglobalvarfilename '4staticoptimization'],'MODELINFO')
    end
    save([OCGRADSOL.basicresultfilename '4staticoptimization'],'sout','gradout')
catch
    failed=1;
end

%-----------------------------------------------------------------
function pathname=datapath()
global OCGRADSOL

pathname=OCGRADSOL.datapath();

% helper function for directional derivative
function dX=directionalderivative(coeff,extremal)
global OCGRADCONT

coeff1=coeff;
coeff1=coeff1+OCGRADCONT.OPTIONS.directionalderivativestep;

[freepar,modelpar]=drearr(coeff1);

extremal1=OCGRADCONT.gradientsolver(extremal,modelpar);
dX=(extremal1.y-extremal.y)/OCGRADCONT.OPTIONS.directionalderivativestep;
