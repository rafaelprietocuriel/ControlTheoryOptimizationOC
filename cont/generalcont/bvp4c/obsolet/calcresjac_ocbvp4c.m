function J=calcresjac_ocbvp4c(tmesh,y,z,freepar,contval,modelpar,odefile,jacobianfile,bc,bcjacobianfile)

global OCMATCONT
domainddata=OCMATCONT.DOMAINDDATA;

F=OCMATCONT.HE.DDATA.F;
% multi-point BVP support
numarc=OCMATCONT.HE.numarc;
leftarcindex=OCMATCONT.HE.TIMEDDATA.leftarcindex;
rightarcindex=OCMATCONT.HE.TIMEDDATA.rightarcindex;
numparameter=OCMATCONT.HE.numparameter;

nonzeroent=OCMATCONT.HE.numdvariables*(2*OCMATCONT.HE.numarc*domainddata(1).numeq+OCMATCONT.HE.numparameter);
ROW=zeros(nonzeroent,1);
COL=zeros(nonzeroent,1);
J=zeros(nonzeroent,1);
In=eye(domainddata.numode);

gen_c=0; %General count for sparse matrix efficiency
ROWcounter=0;
ROWodecounter=sum(OCMATCONT.HE.numboundarycondition);
COLcounter_bR=0;
COLodecounter=0;
for arc=1:numarc
    COLcounter=COLcounter_bR;
    arcindex=OCMATCONT.HE.arcindex(arc);
    FcnArgs={modelpar,freepar,contval,arc};
    numode=domainddata(arcindex).numode;
    odecoord=domainddata(arcindex).odecoord; % coordinate of (first order) odes
    numbc=OCMATCONT.HE.numboundarycondition(arc);
    %----Start Boundary conditions
    numbccoord=1:numbc;

    [yal yar ybl ybr]=determineswitchstate(tmesh,y,z,arc);
    % diffmeshatswitch ... the grid size adjacent to the switching times
    [DRal DRar DRbl DRbr DRp]=bcjacobianfile(yal,yar,ybl,ybr,FcnArgs{:});
    
    if ~isempty(DRal)
        gen_c_start=gen_c+1;
        gen_c=gen_c+numbc*numode;
        J(gen_c_start:gen_c)=DRal(:);

        repmat(ROWcounter+numbccoord.',1,numode);
        ROW(gen_c_start:gen_c)=ans(:);
        repmat(COLcounter+odecoord,numbc,1);
        COL(gen_c_start:gen_c)=ans(:);
        %COLcounter=COLcounter+numode;
    end
    gen_c_start=gen_c+1;
    gen_c=gen_c+numbc*numode;
    J(gen_c_start:gen_c)=DRar(:);
    repmat(ROWcounter+numbccoord.',1,numode);
    ROW(gen_c_start:gen_c)=ans(:);
    repmat(COLcounter+odecoord,numbc,1);
    COL(gen_c_start:gen_c)=ans(:);
    COLcounter=COLcounter+(rightarcindex(arc)-leftarcindex(arc))*numode;

    gen_c_start=gen_c+1;
    gen_c=gen_c+numbc*numode;
    J(gen_c_start:gen_c)=DRbl(:);
    repmat(ROWcounter+numbccoord.',1,numode);
    ROW(gen_c_start:gen_c)=ans(:);
    repmat(COLcounter+odecoord,numbc,1);
    COL(gen_c_start:gen_c)=ans(:);
    COLcounter=COLcounter+numode;
    COLcounter_bR=COLcounter;
    if ~isempty(DRbr)
        gen_c_start=gen_c+1;
        gen_c=gen_c+numbc*numode;
        J(gen_c_start:gen_c)=DRbr(:);
        repmat(ROWcounter+numbccoord.',1,numode);
        ROW(gen_c_start:gen_c)=ans(:);
        repmat(COLcounter+odecoord,numbc,1);
        COL(gen_c_start:gen_c)=ans(:);
    end
    gen_c_start=gen_c+1;
    gen_c=gen_c+numparameter*numbc;
    J(gen_c_start:gen_c)=DRp(:);
    ROW(gen_c_start:gen_c)=repmat(ROWcounter+numbccoord.',numparameter,1);
    repmat(OCMATCONT.HE.coeffcoordmp+(1:numparameter),numbc,1);
    COL(gen_c_start:gen_c)=ans(:);
    ROWcounter=ROWcounter+numbc;
end
for arc=1:numarc
    FcnArgs={modelpar,freepar,contval,arc};

    
    xidx = leftarcindex(arc):rightarcindex(arc);
    xreg = tmesh(xidx);
    yreg = y(:,xidx);
    Freg = F(:,xidx);
    hreg = diff(xreg);

    iidx = xidx(1:end-1);    % mesh interval index
    Nint = length(iidx);
    
    % Collocation equations
    [Ji,dFdpar_i]=jacobianfile(xreg(1),yreg(:,1),FcnArgs{:});

    for i = 1:Nint
        hi = hreg(i);
        % the left mesh point
        xi = xreg(i);
        yi = yreg(:,i);
        Fi = Freg(:,i);
        % the right mesh point
        xip1 = xreg(i+1);
        yip1 = yreg(:,i+1);
        Fip1 = Freg(:,i+1);
        [Jip1, dFdpar_ip1]=jacobianfile(xip1,yip1,FcnArgs{:});

        %the interior points
        xip05 = (xi + xip1)/2;
        yip05 = (yi + yip1)/2 - hi/8*(Fip1 - Fi);

        [Jip05, dFdpar_ip05]=jacobianfile(xip05,yip05,FcnArgs{:});

        twiceJip05 = 2*Jip05;
        % assembly
        gen_c_start=gen_c+1;
        gen_c=gen_c+numode*numode;
        -(In+hi/6*(Ji+twiceJip05*(In+hi/4*Ji)));
        J(gen_c_start:gen_c)=ans(:);
        repmat(ROWodecounter+odecoord.',1,numode);
        ROW(gen_c_start:gen_c)=ans(:);
        repmat(COLodecounter+odecoord,numode,1);
        COL(gen_c_start:gen_c)=ans(:);
        COLodecounter=COLodecounter+numode;
        gen_c_start=gen_c+1;
        gen_c=gen_c+numode*numode;
        In-hi/6*(Jip1+twiceJip05*(In-hi/4*Jip1));
        J(gen_c_start:gen_c)=ans(:);
        repmat(ROWodecounter+odecoord.',1,numode);
        ROW(gen_c_start:gen_c)=ans(:);
        repmat(COLodecounter+odecoord,numode,1);
        COL(gen_c_start:gen_c)=ans(:);

        gen_c_start=gen_c+1;
        gen_c=gen_c+numode*OCMATCONT.HE.numparameter;
        -hi*dFdpar_ip05 + hi^2/12*Jip05* ...
            (dFdpar_ip1-dFdpar_i);
        J(gen_c_start:gen_c)=ans(:);
        repmat(ROWodecounter+odecoord.',1,OCMATCONT.HE.numparameter);
        ROW(gen_c_start:gen_c)=ans(:);
        repmat(OCMATCONT.HE.coeffcoordmp+(1:OCMATCONT.HE.numparameter),numode,1);
        COL(gen_c_start:gen_c)=ans(:);
        ROWodecounter=ROWodecounter+numode;   % next equation
        Ji = Jip1;
        dFdpar_i = dFdpar_ip1;
    end
    COLodecounter=COLodecounter+numode;
end
ROW=ROW(1:gen_c);
COL=COL(1:gen_c);
J=J(1:gen_c);
J=sparse(ROW,COL,J,OCMATCONT.HE.numdvariables-1,OCMATCONT.HE.numdvariables);
