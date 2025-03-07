function [J,ROW,COL]=calcresdynamicjac_bvp4c(tmesh,y,z,freepar,contval,modelpar,jacobianfile,rowoffset)

global OCMATCONT
domainddata=OCMATCONT.DOMAINDDATA;

F=OCMATCONT.HE.DDATA.F;
% multi-point BVP support
numarc=OCMATCONT.HE.numarc;
leftarcindex=OCMATCONT.HE.TIMEDDATA.leftarcindex;
rightarcindex=OCMATCONT.HE.TIMEDDATA.rightarcindex;

nonzeroent=OCMATCONT.HE.numdvariables*(2*domainddata(1).numeq);
ROW=zeros(nonzeroent,1);
COL=zeros(nonzeroent,1);
J=zeros(nonzeroent,1);
In=eye(domainddata.numode);

gen_c=0; %General count for sparse matrix efficiency
ROWcounter=rowoffset;
for arc=1:numarc
    FcnArgs={modelpar,freepar,contval,arc};

    arcindex=OCMATCONT.HE.arcindex(arc);
    numode=domainddata(arcindex).numode;
    odecoord=domainddata(arcindex).odecoord; % coordinate of (first order) odes

    COLcounter=0;
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
        repmat(ROWcounter+odecoord.',1,numode);
        ROW(gen_c_start:gen_c)=ans(:);
        repmat(COLcounter+odecoord,numode,1);
        COL(gen_c_start:gen_c)=ans(:);
        COLcounter=COLcounter+numode;
        gen_c_start=gen_c+1;
        gen_c=gen_c+numode*numode;
        In-hi/6*(Jip1+twiceJip05*(In-hi/4*Jip1));
        J(gen_c_start:gen_c)=ans(:);
        repmat(ROWcounter+odecoord.',1,numode);
        ROW(gen_c_start:gen_c)=ans(:);
        repmat(COLcounter+odecoord,numode,1);
        COL(gen_c_start:gen_c)=ans(:);

        gen_c_start=gen_c+1;
        gen_c=gen_c+numode*OCMATCONT.HE.numparameter;
        -hi*dFdpar_ip05 + hi^2/12*Jip05* ...
            (dFdpar_ip1-dFdpar_i);
        J(gen_c_start:gen_c)=ans(:);
        repmat(ROWcounter+odecoord.',1,OCMATCONT.HE.numparameter);
        ROW(gen_c_start:gen_c)=ans(:);
        repmat(OCMATCONT.HE.coeffcoordmp+(1:OCMATCONT.HE.numparameter),numode,1);
        COL(gen_c_start:gen_c)=ans(:);
        ROWcounter=ROWcounter+numode;   % next equation

        Ji = Jip1;
        dFdpar_i = dFdpar_ip1;
    end
end
ROW=ROW(1:gen_c);
COL=COL(1:gen_c);
J=J(1:gen_c);
% J=sparse(ROW,COL,J,OCMATCONT.HE.numdvariables-1,OCMATCONT.HE.numdvariables);
