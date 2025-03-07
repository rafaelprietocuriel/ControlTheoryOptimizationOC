function out=semismoothoptimization()

out{1}=@generalizedoperatoreq;
out{2}=@generalizedfrechetder;
out{7}=@defaultprocessor;
out{8}=@testfunc;
out{9}=@targetvaluefunc;
out{10}=@probleminit;
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

function res=generalizedoperatoreq(coeff,tangent)
global OCSTATSOL OCSTATCONT
[x,lm,freepar,modelpar]=drearr(coeff);
res=OCSTATSOL.firstordernc(x,lm,modelpar);
res(OCSTATSOL.continuationindex)=0;

function varargout=probleminit(varargin)
coeff=varargin{1};
tangent=varargin{2};
WorkspaceInit(coeff,tangent);

% all done succesfully
varargout{1} = 0;

function J=generalizedfrechetder(coeff,tangent)
global OCSTATSOL OCSTATCONT
[x,lm,freepar,modelpar]=drearr(coeff);
Jx=OCSTATSOL.jacobianfirstordernc(x,lm,modelpar,OCSTATCONT.OPTIONS.generalizedderivativethreshold);
Jp=OCSTATSOL.parameterjacobianfirstordernc(x,lm,modelpar,OCSTATCONT.OPTIONS.generalizedderivativethreshold);

J=[Jx,Jp(:,OCSTATSOL.parametercoordinate)];
% for testing
% numJacOpt.diffvar=1;
% numJacOpt.vectvars=[];
% Jnum=numjaccsd(@generalizedoperatoreq,{coeff,tangent},numel(coeff),numJacOpt);
% Jnum(end,:)=[];
% if max(abs(Jnum(:)-J(:)))>1e-2
%     max(abs(Jnum(:)-J(:)))
% end

%----------------------------------------------------------------
function h=plotcontinuation(coeff,tangent,contdata,makemovie)
global OCSTATSOL OCSTATCONT
[x,lm,freepar,modelpar]=drearr(coeff);
h=OCSTATSOL.plotcontinuation(x,lm,freepar,tangent,modelpar,contdata,makemovie);


%----------------------------------------------------------------
function idx=printcontinuation(coeff,tangent)
global OCSTATCONT OCSTATSOL
idx=[];
if isempty(coeff)
    return
end
if isempty(OCSTATCONT.targettype)
    fprintf(1,' Continuation parameter: %g\n',coeff(OCSTATSOL.continuationindex));
else
    [x,lm,freepar,modelpar]=drearr(coeff);
    switch OCSTATCONT.targettype
        case 'parameter'
            diff=OCSTATCONT.targetvalue-modelpar(OCSTATCONT.targetcoordinate);
        otherwise
            diff=OCSTATCONT.targetvalue-OCSTATSOL.targetfunction(x,lm,modelpar,OCSTATCONT.targettype,OCSTATCONT.targetcoordinate);
    end
    fprintf(1,' Distance to targetvalue: %g\n',diff);
end

%---------------------------------------------------------------
function varargout=defaultprocessor(varargin)
global OCSTATCONT OCSTATSOL

varargout{2}=nan;
% all done succesfully
varargout{1}=0;
%-------------------------------------------------------------
function [out, failed]=testfunc(id,coeff,tangent)
global OCSTATCONT OCSTATSOL

out(1)=0;
failed=[];
%----------------------------------------------------------------
function [out, failed]=targetvaluefunc(id,coeff,tangent)
global OCSTATCONT OCSTATSOL

failed=[];
for ii=id
    switch ii
        case 1
            if isempty(OCSTATCONT.targettype)
                out=[];
            else
                [x,lm,freepar,modelpar]=drearr(coeff);
                switch OCSTATCONT.targettype
                    case 'parameter'
                        out=OCSTATCONT.targetvalue-modelpar(OCSTATCONT.targetcoordinate);
                    otherwise
                        out=OCSTATCONT.targetvalue-OCSTATSOL.targetfunction(x,lm,modelpar,OCSTATCONT.targettype,OCSTATCONT.targetcoordinate);
                end
            end

        otherwise
            out=[];
            failed=1;
    end
end

%----------------------------------------------------------------
function out=formatsolution(coeff,tangent)
global OCSTATCONT OCSTATSOL 
[x,lm,freepar,modelpar]=drearr(coeff);

out.modelparameter=modelpar;
out.parameters=freepar;
out.x=x;
out.lm=lm;
out.modelname=OCSTATSOL.modelname;

%---------------------------------------------------------
function done
WorkspaceDone;
% ---------------------------------------------------------

function WorkspaceInit(coeff,tangent)
global OCSTATCONT OCSTATSOL
% ------------------------------------------------------

function WorkspaceDone

% ------------------------------------------------------

function [x,lm,freepar,modelpar]=drearr(coeff)
global OCSTATCONT OCSTATSOL

modelpar=OCSTATSOL.modelparameter;
freepar=coeff(OCSTATSOL.continuationindex);
x=coeff(OCSTATSOL.variableindex);
lm=coeff(OCSTATSOL.lagrangemulitplierindex);
modelpar(OCSTATSOL.parametercoordinate)=freepar;

%-----------------------------------------------------------------
function failed=saveintermediate(sout,contout,contnum)
global OCSTATCONT OCSTATSOL
failed=0;
MODELINFO.OCSTATCONT=OCSTATCONT;
MODELINFO.OCSTATSOL=OCSTATSOL;
try
    if contnum==1
        save([OCSTATSOL.basicglobalvarfilename '4semismoothoptimization'],'MODELINFO')
    end
    save([OCSTATSOL.basicresultfilename '4semismoothoptimization'],'sout','contout')
catch
    failed=1;
end

%-----------------------------------------------------------------
function pathname=datapath()
global OCSTATSOL

pathname=OCSTATSOL.datapath();
