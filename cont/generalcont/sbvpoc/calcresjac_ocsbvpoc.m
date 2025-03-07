function J=calcresjac_ocsbvpoc(tmesh,y,z,freepar,contval,modelpar,odefile,jacobianfile,bc,bcjacobianfile)
global OCMATCONT
domainddata=OCMATCONT.DOMAINDDATA;

diffmesh=diff(tmesh);
numarc=OCMATCONT.HE.numarc;
leftarcindex=OCMATCONT.HE.TIMEDDATA.leftarcindex;
rightarcindex=OCMATCONT.HE.TIMEDDATA.rightarcindex;

totalnumode=OCMATCONT.HE.totalnumode;
numparameter=OCMATCONT.HE.numparameter;
maxnumcols=max([domainddata.numcols]);
totalnumeq=sum([domainddata(OCMATCONT.HE.arcindex).numode])+numarc-1;

nonzeroent=(totalnumode+numparameter)*(2*totalnumode+maxnumcols*totalnumeq+numparameter)+...
    (rightarcindex(end)-1)*maxnumcols*totalnumode*totalnumeq+...
    (rightarcindex(end)-1)*maxnumcols*totalnumeq*maxnumcols*totalnumeq+...
    (rightarcindex(end)-1)*2*totalnumeq*totalnumeq+...
    (rightarcindex(end)-1)*totalnumode*totalnumeq*maxnumcols+...
    (rightarcindex(end)-1)*(totalnumode+maxnumcols*totalnumeq)*numparameter;
ROW=zeros(nonzeroent,1);
COL=zeros(nonzeroent,1);
J=zeros(nonzeroent,1);
numbc=OCMATCONT.HE.totalnumboundarycondition;
numbccoord=1:numbc;

ROWcounterCols=numbc;
COLcounterCols=0;
[DRa DRb DRp]=bcjacobianfile(OCMATCONT.HE.yL,OCMATCONT.HE.yR,modelpar,freepar,contval);

gen_c=0; %General count for sparse matrix efficiency
gen_c_start=gen_c+1;
gen_c=gen_c+numparameter*numbc;
J(gen_c_start:gen_c)=DRp(:);
ROW(gen_c_start:gen_c)=repmat(numbccoord.',numparameter,1);
repmat(OCMATCONT.HE.coeffcoordmp+(1:numparameter),numbc,1);
COL(gen_c_start:gen_c)=ans(:);

COLcounter=0;
COLBCcounter=0;
for arc=1:numarc
    arcindex=OCMATCONT.HE.arcindex(arc);
    idx=leftarcindex(arc)-arc+1:rightarcindex(arc)-arc;
    collocationpts=OCMATCONT.HE.DDATA(arc).collocationpoints;
    psival=domainddata(arcindex).psival;
    psival0=domainddata(arcindex).psival0;
    numcols=domainddata(arcindex).numcols;
    numcolscoord=domainddata(arcindex).numcolscoord;
    numeq=domainddata(arcindex).numeq;
    numae=domainddata(arcindex).numae;
    numode=domainddata(arcindex).numode;
    odecoord=domainddata(arcindex).odecoord; % coordinate of (first order) odes
    aecoord=domainddata(arcindex).aecoord; % coordinate of algebraic equations
    %diffmeshatswitch=OCMATCONT.HE.gridatswitch(arc,1:4);
    diffmesh_arc=diffmesh(leftarcindex(arc):rightarcindex(arc)-1);
    N=rightarcindex(arc)-leftarcindex(arc)-1;
    %----Start Boundary conditions

    COLBCcounterstart=COLBCcounter+1;
    COLBCcounter=COLBCcounter+numeq;
    DRar=DRa(numbccoord,COLBCcounterstart:COLBCcounter);
    DRbl=DRb(numbccoord,COLBCcounterstart:COLBCcounter);
    gen_c_start=gen_c+1;
    gen_c=gen_c+numbc*numode;
    DRar(numbccoord,odecoord);
    J(gen_c_start:gen_c)=ans(:);
    repmat(numbccoord.',1,numode);
    ROW(gen_c_start:gen_c)=ans(:);
    repmat(COLcounter+odecoord,numbc,1);
    COL(gen_c_start:gen_c)=ans(:);
    COLcounter=COLcounter+(rightarcindex(arc)-leftarcindex(arc)-1)*(numcols*numeq+numode);

    gen_c_start=gen_c+1;
    gen_c=gen_c+numbc*numode;
    DRbl(numbccoord,odecoord);
    J(gen_c_start:gen_c)=ans(:);
    repmat(numbccoord.',1,numode);
    ROW(gen_c_start:gen_c)=ans(:);
    repmat(COLcounter+odecoord,numbc,1);
    COL(gen_c_start:gen_c)=ans(:);
    COLcounter=COLcounter+numode;
    for ll=numcolscoord
        gen_c_start=gen_c+1;
        gen_c=gen_c+numbc*numode;
        diffmesh_arc(N)*DRbl(numbccoord,odecoord)*psival(1,ll,numcols+1);
        J(gen_c_start:gen_c)=ans(:);
        repmat(numbccoord.',1,numode);
        ROW(gen_c_start:gen_c)=ans(:);
        repmat(COLcounter+odecoord,numbc,1);
        COL(gen_c_start:gen_c)=ans(:);
        COLcounter=COLcounter+numode;
        if ~isempty(aecoord)
            gen_c_start=gen_c+1;
            gen_c=gen_c+numbc*numae;
            DRbl(numbccoord,aecoord)*psival0(1,ll,numcols+1);
            J(gen_c_start:gen_c)=ans(:);
            ROW(gen_c_start:gen_c)=repmat(numbccoord.',numae,1);
            repmat(COLcounter+(1:numae),numbc,1);
            COL(gen_c_start:gen_c)=ans(:);
            COLcounter=COLcounter+numae;
        end
    end
    
    
    %----------Start Jacobian entries at collocations and mesh points
    eqcoord=domainddata(arcindex).eqcoord;
    eqderivative=eye(numode);
    if numae
        eqderivative=blkdiag(eye(numode),zeros(numae));
    end
    fz=zeros(numeq,2,numcols);
    for step=leftarcindex(arc):rightarcindex(arc)-1
        countermeshpts=step-leftarcindex(arc)+1;
        for kk=numcolscoord
            fz(odecoord,1,numcolscoord)=evaluatePcols(countermeshpts,diffmesh_arc,y(odecoord,1,idx),z(odecoord,numcolscoord,idx),odecoord,psival,numcols);
            fz(odecoord,2,numcolscoord)=z(odecoord,numcolscoord,idx(countermeshpts));
            fz(aecoord,1,numcolscoord)=z(aecoord,numcolscoord,idx(countermeshpts));
            %ii indicates, which component of g_ik will be differentiated with respect
            %to y und z
            [Dg Dpg]=jacobianfile(collocationpts(countermeshpts,kk),fz(:,1,kk),modelpar,freepar,contval,arc);
            for ii=eqcoord
                ROWcounterCols=ROWcounterCols+1;
                %Differentiation with respect to y_i, sequence: y_i1 to
                %y_i(maxorder)
                counter=0;
                gen_c_start=gen_c+1;
                gen_c=gen_c+numode;
                J(gen_c_start:gen_c)=-Dg(ii,odecoord);
                ROW(gen_c_start:gen_c)=ROWcounterCols;
                COL(gen_c_start:gen_c)=COLcounterCols+odecoord;
                counter=counter+numode+1;
                help=zeros(numeq,numcols);
                for jj=odecoord
                    help(jj,numcolscoord)=eqderivative(ii,jj).*(kk==(numcolscoord))-Dg(ii,jj)*diffmesh_arc(countermeshpts).*psival(1,numcolscoord,kk);
                end
                for jj=aecoord
                    help(jj,numcolscoord)=eqderivative(ii,jj).*(kk==(numcolscoord))-Dg(ii,jj).*(kk==(numcolscoord));
                end
                col=COLcounterCols+counter:COLcounterCols+counter+numeq*numcols-1;
                gen_c_start=gen_c+1;
                gen_c=gen_c+numeq*numcols;
                ROW(gen_c_start:gen_c)=ROWcounterCols;
                COL(gen_c_start:gen_c)=col;
                J(gen_c_start:gen_c)=help(:).';
                if numparameter>0
                    gen_c_start=gen_c+1;
                    gen_c=gen_c+numparameter;
                    ROW(gen_c_start:gen_c)=ROWcounterCols;
                    COL(gen_c_start:gen_c)=OCMATCONT.HE.coeffcoordmp+[1:numparameter];
                    J(gen_c_start:gen_c)=-Dpg(ii,1:numparameter);
                end
            end
        end
        if step<rightarcindex(arc)-1
            gen_c_start=gen_c+1;
            gen_c=gen_c+numode;
            ROW(gen_c_start:gen_c)=ROWcounterCols+odecoord;
            COL(gen_c_start:gen_c)=COLcounterCols+odecoord;
            J(gen_c_start:gen_c)=1;
            gen_c_start=gen_c+1;
            gen_c=gen_c+numode;
            ROW(gen_c_start:gen_c)=ROWcounterCols+odecoord;
            COL(gen_c_start:gen_c)=COLcounterCols+numcols*numeq+numode+odecoord;
            J(gen_c_start:gen_c)=-1;
            gen_c_start=gen_c+1;
            gen_c=gen_c+numcols*numode;
            ROWcounterCols+repmat(odecoord,numcols,1);
            ROW(gen_c_start:gen_c)=ans(:).';
            COLcounterCols+numode+repmat([0:numcols-1]'*numeq,1,numode)+repmat(odecoord,numcols,1);
            COL(gen_c_start:gen_c)=ans(:);
            repmat(diffmesh_arc(countermeshpts)*psival(1,numcolscoord,numcols+1).',1,numode);
            J(gen_c_start:gen_c)=ans(:);
            ROWcounterCols=ROWcounterCols+numode;
        end
        COLcounterCols=COLcounterCols+numcols*numeq+numode;
    end
end
ROW=ROW(1:gen_c);
COL=COL(1:gen_c);
J=J(1:gen_c);
J=sparse(ROW,COL,J,OCMATCONT.HE.numdvariables-1,OCMATCONT.HE.numdvariables);
