function out=indifferencegradsolution()

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
out{26}=@saveintermediate;
out{27}=@datapath;
out{30}=@printcontinuation;

function [res,extremal,totgraditer,totdgraditer,totlineiter,totdlineiter,max_numrepdd,statinfo]=operatoreq(coeff,extremal,tangent)
global OCGRADCONT OCGRADSOL
[freepar,modelpar]=drearr(coeff);
res=zeros(OCGRADCONT.HE.numdvariables,1);
y0=OCGRADCONT.initialstatevalue;
switch OCGRADCONT.conttype
    case 'initialstate'
        if OCGRADCONT.numfreestatevector==1
            y0=y0+freepar(OCGRADCONT.continuationcoordinate)*OCGRADCONT.freestatevector;
        elseif OCGRADCONT.numfreestatevector==2
            y0=y0+freepar(OCGRADCONT.continuationcoordinate)*OCGRADCONT.freestatevector(:,2)+freepar(OCGRADCONT.freestatevectorcoordinate)*OCGRADCONT.freestatevector(:,1);
        end
    case 'parameter'
        if OCGRADCONT.numfreestatevector==1
            y0=y0+freepar(OCGRADCONT.freestatevectorcoordinate)*OCGRADCONT.freestatevector;
        end
end
totgraditer=0;
totlineiter=0;
for ii=1:OCGRADSOL.degree
    extremal(ii).y=y0;
    if strcmp(OCGRADCONT.trajectoryclass{ii},'inf')
        OCGRADCONT.infinitetimeendconditions=1;
        extremal(ii).cst_y(:,end)=freepar(OCGRADCONT.freecostatecoordinate{ii});
    else
        OCGRADCONT.infinitetimeendconditions=0;
    end
    [extremal(ii),graditer,lineiter,tmpstatinfo]=OCGRADCONT.gradientsolver(extremal(ii),modelpar);
    if ii==1
        statinfo=tmpstatinfo;
    else
        if ~isempty(statinfo)
            statinfo=[statinfo tmpstatinfo];
        end
    end
    totgraditer=totgraditer+graditer;
    totlineiter=totlineiter+lineiter;
end

totdgraditer=0;
totdlineiter=0;
max_numrepdd=-inf;
rescounter=0;
for ii=1:OCGRADSOL.degree % last parameter is the continuaton parameter
    if strcmp(OCGRADCONT.trajectoryclass{ii},'inf')
        OCGRADCONT.infinitetimeendconditions=1;
    else
        OCGRADCONT.infinitetimeendconditions=0;
    end
    [extremal(ii).tangent,dgraditer,dlineiter,numrepdd]=directionalderivative(coeff,extremal(ii));
    totdgraditer=totdgraditer+dgraditer;
    totdlineiter=totdlineiter+dlineiter;
    max_numrepdd=max(max_numrepdd,numrepdd);
    if strcmp(OCGRADCONT.trajectoryclass{ii},'inf')
        res(rescounter+(1:OCGRADCONT.freecostatenum{ii}))=OCGRADSOL.asymptoticmatrix{ii}.'*(OCGRADSOL.saddlepoint{ii}-[extremal(ii).y(:,OCGRADCONT.TIMEMESH.num);freepar(OCGRADCONT.freecostatecoordinate{ii})]);
        rescounter=rescounter+OCGRADCONT.freecostatenum{ii};
    end
    if ii>1
        if strcmp(OCGRADCONT.trajectoryclass{ii},'inf')
            res(rescounter+1)=OCGRADSOL.asymptoticobjectivevalue(0,extremal(ii).y(:,1),extremal(ii).v(:,1),extremal(ii).cst_y(:,1),modelpar)- ...
                OCGRADSOL.asymptoticobjectivevalue(0,extremal(ii-1).y(:,1),extremal(ii-1).v(:,1),extremal(ii-1).cst_y(:,1),modelpar);
        else
            res(rescounter+1)=extremal(ii).objectivevalue-extremal(ii-1).objectivevalue;
        end
        rescounter=rescounter+1;
    end
end
res(OCGRADCONT.HE.numdvariables,1)=0;

function res=operatoreqred(coeff,extremal,tangent)
global OCGRADCONT OCGRADSOL
[freepar,modelpar]=drearr(coeff);
res=zeros(OCGRADCONT.HE.numdvariables-1,1);
y0=OCGRADCONT.initialstatevalue;
switch OCGRADCONT.conttype
    case 'initialstate'
        if OCGRADCONT.numfreestatevector==1
            y0=y0+freepar(OCGRADCONT.continuationcoordinate)*OCGRADCONT.freestatevector;
        elseif OCGRADCONT.numfreestatevector==2
            y0=y0+freepar(OCGRADCONT.continuationcoordinate)*OCGRADCONT.freestatevector(:,2)+freepar(OCGRADCONT.freestatevectorcoordinate)*OCGRADCONT.freestatevector(:,1);
        end
    case 'parameter'
        if OCGRADCONT.numfreestatevector==1
            y0=y0+freepar(OCGRADCONT.freestatevectorcoordinate)*OCGRADCONT.freestatevector;
        end
end
rescounter=0;
for ii=1:OCGRADSOL.degree
    extremal(ii).y=y0;
    if strcmp(OCGRADCONT.trajectoryclass{ii},'inf')
        OCGRADCONT.infinitetimeendconditions=1;
    else
        OCGRADCONT.infinitetimeendconditions=0;
    end
    extremal(ii)=OCGRADCONT.gradientsolver(extremal(ii),modelpar);
    if strcmp(OCGRADCONT.trajectoryclass{ii},'inf')
        res(rescounter+(1:OCGRADCONT.freecostatenum{ii}))=OCGRADSOL.asymptoticmatrix{ii}.'*(OCGRADSOL.saddlepoint{ii}-[extremal(ii).y(:,OCGRADCONT.TIMEMESH.num);freepar(OCGRADCONT.freecostatecoordinate{ii})]);
        rescounter=rescounter+OCGRADCONT.freecostatenum{ii};
    end
    if ii>1
        if strcmp(OCGRADCONT.trajectoryclass{ii},'inf')
            res(rescounter+1)=OCGRADSOL.asymptoticobjectivevalue(0,extremal(ii).y(:,1),extremal(ii).v(:,1),extremal(ii).cst_y(:,1),modelpar)- ...
                OCGRADSOL.asymptoticobjectivevalue(0,extremal(ii-1).y(:,1),extremal(ii-1).v(:,1),extremal(ii-1).cst_y(:,1),modelpar);
        else
            res(rescounter+1)=extremal(ii).objectivevalue-extremal(ii-1).objectivevalue;
        end
        rescounter=rescounter+1;
    end
end

function J=frechetder(coeff,extremal,tangent,varargin)
global OCGRADCONT OCGRADSOL
% [freepar,modelpar]=drearr(coeff);
% % derivative d/dp V = int_0^Td/dp g(y(p),v(p))dt + d/dpS
% for ii=1:OCGRADSOL.degree
%     dSdy=OCGRADSOL.DsalvagevalueDy(extremal(ii).timehorizon,extremal(ii).y(:,end),modelpar);
%     dSdp=OCGRADSOL.DsalvagevalueDp(extremal(ii).timehorizon,extremal(ii).y(:,end),modelpar);
%     dO(:,ii)=extremal(ii).tangent.o(:,end);
%     dy_end=extremal(ii).tangent.y(end-OCGRADCONT.state_num.concentrated+1:end,end);
%     if ~isempty(OCGRADCONT.freeparameterindex)
%         dO(OCGRADCONT.freeparametercoordinate,ii)=dO(OCGRADCONT.freeparametercoordinate,ii)+ ...
%             dSdy.'*dy_end+dSdp(OCGRADCONT.freeparameterindex);
%     end
% end
% 
% J0=dO(:,2).'-dO(:,1).';

%%%%
J=zeros(OCGRADCONT.HE.numdvariables-1,OCGRADCONT.HE.numdvariables);
res0=operatoreqred(coeff,extremal,tangent);
directionalderivativestep=OCGRADCONT.OPTIONS.directionalderivativestep;
for ii=1:length(coeff) % loop for the different solutions
    coeff1=coeff;
    coeff1(ii)=coeff1(ii)+directionalderivativestep;
    res=operatoreqred(coeff1,extremal,tangent);
    J(:,ii)=(res-res0)/directionalderivativestep;
end

function varargout=probleminit(varargin)
coeff=varargin{1};
tangent=varargin{2};
WorkspaceInit(coeff,tangent);

% all done succesfully
varargout{1} = 0;

%-------------------------------------------------------------------------
function [coeff,extremal,tangent,graditer]=gradientsolution(coeff,extremal,tangent)
coeff=[];
extremal=[];
tangent=[];
graditer=[];

%-------------------------------------------------------------------------
function [precoeff,extremal]=predictextremal(coeff,extremal,stepwidth,tangent)
global OCGRADCONT OCGRADSOL
tot_extremal_tangent=[];
tot_extremal_coeff=[];
for ii=1:OCGRADSOL.degree
    tot_extremal_tangent=[tot_extremal_tangent;extremal(ii).tangent.coeff(:,end)];
    tot_extremal_coeff=[tot_extremal_coeff;extremal(ii).coeff];
end
tot_extremal_tangent=tot_extremal_tangent*tangent(end);
tot_extremal_tangent=[tot_extremal_tangent;tangent];
tot_extremal_tangent=tot_extremal_tangent/norm(tot_extremal_tangent);
tot_extremal_coeff=[tot_extremal_coeff;coeff];
tot_precoeff=tot_extremal_coeff+stepwidth*tot_extremal_tangent;
precoeff=tot_precoeff((end-OCGRADCONT.HE.numdvariables+1):end);
for ii=1:OCGRADSOL.degree
    extremal(ii).y=tot_precoeff(OCGRADCONT.initiallocstatecoordinate{ii});
    extremal(ii).cst_y=tot_precoeff(OCGRADCONT.endloccostatecoordinate{ii});
    if OCGRADCONT.control_num.concentrated==1
        extremal(ii).v=tot_precoeff(OCGRADCONT.loccontrolcoordinate{ii}).';
    else
        extremal(ii).v=tot_precoeff(OCGRADCONT.loccontrolcoordinate{ii});
    end
    extremal(ii).tangent=[];
end


%-------------------------------------------------------------------------
function [coeffpre,extremal]=predictextremaldiff(coeff1,coeff2,extremal1,extremal2,fac)
global OCGRADCONT OCGRADSOL
coeffpre=coeff1+fac*(coeff2-coeff1);
extremal=extremal1;

for ii=1:OCGRADSOL.degree
    extremal(ii).v=extremal1(ii).v+fac*(extremal2(ii).v-extremal1(ii).v);
    extremal(ii).y=extremal1(ii).y+fac*(extremal2(ii).y-extremal1(ii).y);
    extremal(ii).cst_y=extremal1(ii).cst_y+fac*(extremal2(ii).cst_y-extremal1(ii).cst_y);
    extremal(ii).tangent=[];
    extremal(ii).dv=[];
end
%----------------------------------------------------------------
function h=plotcontinuation(coeff,extremal,tangent,contdata,makemovie)
global OCGRADSOL OCGRADCONT
[freepar,modelpar]=drearr(coeff);
h=OCGRADSOL.plotcontinuation(freepar,extremal,tangent,modelpar,contdata,OCGRADCONT.problem_func,makemovie);


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
    switch OCGRADCONT.conttype
        case 'initialstate'
            fprintf(1,' distance from target value: %g\n',OCGRADCONT.targetvalue-extremal(1).y(OCGRADCONT.targetindex,1));
        case 'parameter'
            freepar=drearr(coeff);
            fprintf(1,' distance from target value: %g\n',OCGRADCONT.targetvalue-freepar(OCGRADCONT.continuationcoordinate));
    end
    
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
                    case 'initialstate'
                        out=OCGRADCONT.targetvalue-extremal(1).y(OCGRADCONT.targetindex,1);
                    case 'parameter'
                        freepar=drearr(coeff);
                        out=OCGRADCONT.targetvalue-freepar(OCGRADCONT.continuationcoordinate);
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
b=zeros(1,OCGRADSOL.degree);
for ii=1:OCGRADSOL.degree
    val=OCGRADSOL.admissible(extremal(ii).t,extremal(ii).y,extremal(ii).v,modelpar);
    b(ii)=any(val(:)<0);
end
b=any(b);

%----------------------------------------------------------------
function out=formatsolution(coeff,extremal,tangent)
global OCGRADCONT OCGRADSOL
[freepar,modelpar]=drearr(coeff);

for ii=1:OCGRADSOL.degree
    extremal(ii).modelname=OCGRADCONT.modelname;
    extremal(ii).modelparameter=modelpar;
    extremal(ii).solver.gradient=OCGRADCONT.gradientsolver;
    extremal(ii).solver.ode=OCGRADCONT.gradientodesolver;
    extremal(ii).solver.option.R0=OCGRADCONT.OPTIONS.R0;
    extremal(ii).solver.option.gradientmappingmethod=OCGRADCONT.OPTIONS.gradientmappingmethod;
    extremal(ii).solver.option.gradtol=OCGRADCONT.OPTIONS.gradtol;
    extremal(ii).solver.option.gamma=OCGRADCONT.OPTIONS.gamma;
    extremal(ii).solver.option.mu=OCGRADCONT.OPTIONS.mu;
    extremal(ii).contclass=OCGRADCONT.problem_func;
    extremal(ii).conttype=OCGRADCONT.conttype;
    extremal(ii).contindex=OCGRADCONT.contindex;
    extremal(ii).freeparameterindex=OCGRADCONT.freeparameterindex;
    extremal(ii).freestatevector=OCGRADCONT.freestatevector;
    extremal(ii).trajectoryclass=OCGRADCONT.trajectoryclass;
    extremal(ii).freeparametercoordinate=OCGRADCONT.freeparametercoordinate;
    if strcmp(OCGRADCONT.trajectoryclass{ii},'inf')
        extremal(ii).freestatevectorcoordinate=OCGRADCONT.freestatevectorcoordinate;
    end
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
freepar=coeff;
switch OCGRADCONT.conttype
    case 'initialstate'
        if ~isempty(OCGRADCONT.freeparameterindex)
            modelpar(OCGRADCONT.freeparameterindex)=coeff(OCGRADCONT.freeparametercoordinate);
        end
    case 'parameter'
        if ~isempty(OCGRADCONT.freeparameterindex)
            modelpar([OCGRADCONT.freeparameterindex OCGRADCONT.contindex])=coeff([OCGRADCONT.freeparametercoordinate OCGRADCONT.continuationcoordinate]);
        else
            modelpar(OCGRADCONT.contindex)=coeff(OCGRADCONT.continuationcoordinate);
        end
end

%-----------------------------------------------------------------
function failed=saveintermediate(sout,gradout,contnum)
global OCGRADCONT OCGRADSOL
failed=0;
MODELINFO.OCGRADCONT=OCGRADCONT;
MODELINFO.OCGRADSOL=OCGRADSOL;
try
    if contnum==1
        save([OCGRADSOL.basicglobalvarfilename '4indifferencegradsolution'],'MODELINFO')
    end
    save([OCGRADSOL.basicresultfilename '4indifferencegradsolution'],'sout','gradout')
catch
    failed=1;
end

%-----------------------------------------------------------------
function pathname=datapath()
global OCGRADSOL

pathname=OCGRADSOL.datapath();

% helper function for directional derivative
function [dX,totgraditer,totlineiter,ctr_max]=directionalderivative(coeff,extremal)
global OCGRADSOL OCGRADCONT

%dX=OCGRADCONT.initdX;
ctrder=0;
y=extremal.y;
v=extremal.v;
cst_y=extremal.cst_y;
o=extremal.o;
ctr_max=-inf;
dX.y=[];
dX.v=[];
dX.cst_y=[];
dX.o=[];
dX.coeff=[];
totgraditer=0;
totlineiter=0;
for ii=1:length(coeff) % loop for the different solutions
    ctrder=ctrder+1;
    directionalderivativestep=OCGRADCONT.OPTIONS.directionalderivativestep;
    ctr=0;
    while 1
        ctr=ctr+1;
        coeff1=coeff;
        coeff1(ii)=coeff1(ii)+directionalderivativestep;
        [freepar,modelpar]=drearr(coeff1);
        y0=OCGRADCONT.initialstatevalue;
        switch OCGRADCONT.conttype
            case 'initialstate'
            case 'initialstate'
                if OCGRADCONT.numfreestatevector==1
                    y0=y0+freepar(OCGRADCONT.continuationcoordinate)*OCGRADCONT.freestatevector;
                elseif OCGRADCONT.numfreestatevector==2
                    y0=y0+freepar(OCGRADCONT.continuationcoordinate)*OCGRADCONT.freestatevector(:,2)+freepar(OCGRADCONT.freestatevectorcoordinate)*OCGRADCONT.freestatevector(:,1);
                end
            case 'parameter'
                if OCGRADCONT.numfreestatevector==1
                    y0=y0+freepar(OCGRADCONT.freestatevectorcoordinate)*OCGRADCONT.freestatevector;
                end
        end
        extremal.y=y0;
        [extremal1,graditer,lineiter]=OCGRADCONT.gradientsolver(extremal,modelpar);
        totgraditer=totgraditer+graditer;
        totlineiter=totlineiter+lineiter;
        dx.y=(extremal1.y-y)/directionalderivativestep;
        dx.v=(extremal1.v-v)/directionalderivativestep;
        dx.cst_y=(extremal1.cst_y-cst_y)/directionalderivativestep;
        dx.o=(extremal1.o-o)/directionalderivativestep;
        dx.coeff=[dx.v(:);dx.y(:,1);dx.cst_y(:,OCGRADCONT.TIMEMESH.num)];
        if norm(dx.coeff)~=0
            ctr_max=max(ctr_max,ctr);
            dX.y=[dX.y;dx.y];
            dX.cst_y=[dX.cst_y;dx.cst_y];
            dX.v=[dX.v;dx.v];
            dX.o=[dX.o;dx.o];
            dX.coeff=[dX.coeff,dx.coeff];
            break
        end
        if ctr>0
            ctr_max=max(ctr_max,ctr);
            dX.y=[dX.y;dx.y];
            dX.cst_y=[dX.cst_y;dx.cst_y];
            dX.v=[dX.v;dx.v];
            dX.o=[dX.o;dx.o];
            dX.coeff=[dX.coeff,dx.coeff];
            break
        end
        directionalderivativestep=directionalderivativestep*5;
    end
end

