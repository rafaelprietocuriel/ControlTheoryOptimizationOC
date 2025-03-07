function [J,ROW,COL,ROWcounter]=calcresbcjac_bvp4c(tmesh,y,z,freepar,contval,modelpar,bcjacobianfile)

global OCMATCONT
domainddata=OCMATCONT.DOMAINDDATA;

% multi-point BVP support
numarc=OCMATCONT.HE.numarc;
leftarcindex=OCMATCONT.HE.TIMEDDATA.leftarcindex;
rightarcindex=OCMATCONT.HE.TIMEDDATA.rightarcindex;
numparameter=OCMATCONT.HE.numparameter;

nonzeroent=OCMATCONT.HE.numdvariables*OCMATCONT.HE.numparameter;
ROW=zeros(nonzeroent,1);
COL=zeros(nonzeroent,1);
J=zeros(nonzeroent,1);

gen_c=0; %General count for sparse matrix efficiency
ROWcounter=0;
COLcounter=0;
for arc=1:numarc
    FcnArgs={modelpar,freepar,contval,arc};

    arcindex=OCMATCONT.HE.arcindex(arc);
    numode=domainddata(arcindex).numode;
    numparametermc=OCMATCONT.HE.numparametermc;
    odecoord=domainddata(arcindex).odecoord; % coordinate of (first order) odes

    %----Start Boundary conditions
    if arc==1
        numbc=OCMATCONT.HE.numinitialcondition;
    else
        numbc=0;
    end
    if numarc>1
        if arc>1
            numbc=numbc+domainddata(OCMATCONT.HE.arcindex(arc-1)).numode;
        else
        end
        if arc<numarc
            numbc=numbc+1;
        else
            numbc=numbc+OCMATCONT.HE.numendcondition;
        end
    else
        numbc=numbc+OCMATCONT.HE.numendcondition;
    end
    numbccoord=1:numbc;

    [yal yar ybl ybr]=determineswitchstate(tmesh,y,z,arc);
    % diffmeshatswitch ... the grid size adjacent to the switching times
    [DRal DRar DRbl DRbr DRp]=bcjacobianfile(yal,yar,ybl,ybr,FcnArgs{:});
    
    if ~isempty(DRal)
        gen_c_start=gen_c+1;
        gen_c=gen_c+numbc*(numbc+numparametermc);
        J(gen_c_start:gen_c)=DRal(:);

        repmat(ROWcounter+(1:numbc+numparametermc)',1,numbc);
        ROW(gen_c_start:gen_c)=ans(:);
        repmat(COLcounter+odecoord,numbc,1);
        COL(gen_c_start:gen_c)=ans(:);
        COLcounter=COLcounter+numode;
    end
    gen_c_start=gen_c+1;
    gen_c=gen_c+numbc*(numbc+numparametermc);
    J(gen_c_start:gen_c)=DRar(:);
    repmat(ROWcounter+(1:numbc+numparametermc)',1,numbc);
    ROW(gen_c_start:gen_c)=ans(:);
    repmat(COLcounter+odecoord,numbc+numparametermc,1);
    COL(gen_c_start:gen_c)=ans(:);
    COLcounter=COLcounter+(rightarcindex(arc)-leftarcindex(arc))*numode;

    gen_c_start=gen_c+1;
    gen_c=gen_c+numbc*(numbc+numparametermc);
    J(gen_c_start:gen_c)=DRbl(:);
    repmat(ROWcounter+(1:numbc+numparametermc)',1,numbc);
    ROW(gen_c_start:gen_c)=ans(:);
    repmat(COLcounter+odecoord,numbc+numparametermc,1);
    COL(gen_c_start:gen_c)=ans(:);
    COLcounter=COLcounter+numode;
    if ~isempty(DRbr)
        gen_c_start=gen_c+1;
        gen_c=gen_c+numbc*(numbc+numparametermc);
        J(gen_c_start:gen_c)=DRbr(:);
        repmat(ROWcounter+(1:numbc+numparametermc)',1,numbc);
        ROW(gen_c_start:gen_c)=ans(:);
        repmat(COLcounter+odecoord,numbc+numparametermc,1);
        COL(gen_c_start:gen_c)=ans(:);
    end
    gen_c_start=gen_c+1;
    gen_c=gen_c+numparameter*(numbc+numparametermc);
    J(gen_c_start:gen_c)=DRp(:);
    repmat(ROWcounter+(1:numbc+numparametermc)',1,numparameter);
    ROW(gen_c_start:gen_c)=ans(:);
    repmat(OCMATCONT.HE.coeffcoordmp+(1:numparameter),numbc+numparametermc,1);
    COL(gen_c_start:gen_c)=ans(:);
    ROWcounter=ROWcounter+numbc+numparametermc;
end
ROW=ROW(1:gen_c);
COL=COL(1:gen_c);
J=J(1:gen_c);