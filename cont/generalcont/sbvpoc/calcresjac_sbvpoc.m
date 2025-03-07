function J=calcresjac_sbvpoc(tmesh,y,z,freepar,modelpar,ode,bc,odejac,bcjac)
global OCMATCONT OCBVP
domainddata=OCMATCONT.DOMAINDDATA;

diffmesh=diff(tmesh);
numarc=OCMATCONT.HE.numarc;
leftarcindex=OCMATCONT.HE.TIMEDDATA.leftarcindex;
rightarcindex=OCMATCONT.HE.TIMEDDATA.rightarcindex;
Lidx = OCBVP.Lidx;
Ridx = OCBVP.Ridx;
FcnArgs = {0,freepar,modelpar};    % Pass the region index to the ODE function.

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

ROWcounterCols=totalnumeq;
COLcounterCols=0;


if isempty(bcjac)   % use numerical approx
    [dGdya,dGdyb,nbc,dGdpar] = BCnumjaccol(bc,tmesh,y(:,1,Lidx),y(:,1,Ridx),z(:,:,Lidx),z(:,:,Ridx),FcnArgs(2:3));
else  % use analytical Jacobian
    [dGdya,dGdyb,dGdpar] = bcjac(tmesh,y(:,1,Lidx),y(:,1,Ridx),z(:,:,Lidx),z(:,:,Ridx),FcnArgs{2:3});
end
last_cols(1:nBCs,:) = dGdpar;

gen_c=0; %General count for sparse matrix efficiency
ROWcounter=0;
COLcounter=0;
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

    %----Start Boundary conditions
%     if arc==1
%         numbc=OCMATCONT.HE.numinitialcondition;
%     else
%         numbc=0;
%     end
    numbc=OCMATCONT.HE.numboundarycondition(arc);
    if numarc>1
        if arc>1
            psivalL=domainddata(OCMATCONT.HE.arcindex(arc-1)).psival;
            psival0L=domainddata(OCMATCONT.HE.arcindex(arc-1)).psival0;
            numcolsL=domainddata(OCMATCONT.HE.arcindex(arc-1)).numcols;
            numcolscoordL=domainddata(OCMATCONT.HE.arcindex(arc-1)).numcolscoord;
            numaeL=domainddata(OCMATCONT.HE.arcindex(arc-1)).numae;
            numodeL=domainddata(OCMATCONT.HE.arcindex(arc-1)).numode;
            odecoordL=domainddata(OCMATCONT.HE.arcindex(arc-1)).odecoord; % coordinate of (first order) odes
            aecoordL=domainddata(OCMATCONT.HE.arcindex(arc-1)).aecoord; % coordinate of (first order) odes
        end
        if arc<numarc
            psivalR=domainddata(OCMATCONT.HE.arcindex(arc+1)).psival;
            psival0R=domainddata(OCMATCONT.HE.arcindex(arc+1)).psival0;
            numcolsR=domainddata(OCMATCONT.HE.arcindex(arc+1)).numcols;
            numcolscoordR=domainddata(OCMATCONT.HE.arcindex(arc+1)).numcolscoord;
            numaeR=domainddata(OCMATCONT.HE.arcindex(arc+1)).numae;
            numodeR=domainddata(OCMATCONT.HE.arcindex(arc+1)).numode;
            odecoordR=domainddata(OCMATCONT.HE.arcindex(arc+1)).odecoord; % coordinate of (first order) odes
            aecoordR=domainddata(OCMATCONT.HE.arcindex(arc+1)).aecoord; % coordinate of (first order) odes
        end
    end
    numbccoord=1:numbc;

    [yal yar ybl ybr diffmeshatswitch]=determineswitchstate(tmesh,y,z,arc);
    % diffmeshatswitch ... the grid size adjacent to the switching times
    [DRal DRar DRbl DRbr DRp]=bcjac(yal,yar,ybl,ybr,modelpar,freepar,contval,arc);

    if ~isempty(DRal)
        gen_c_start=gen_c+1;
        gen_c=gen_c+numbc*numodeL;
        % derivative with respect to y
        DRal(numbccoord,odecoordL);
        J(gen_c_start:gen_c)=ans(:);
        
        ROW(gen_c_start:gen_c)=repmat(ROWcounter+numbccoord.',numodeL,1);
        repmat(COLcounter+odecoordL,numbc,1);
        COL(gen_c_start:gen_c)=ans(:);
        COLcounter=COLcounter+numodeL;
        for ll=numcolscoordL
            gen_c_start=gen_c+1;
            gen_c=gen_c+numbc*numodeL;
            %diffmeshatswitch(1)*DRal(numbccoord,odecoordL).*psivalL(ones(numbc,1),ll(1,ones(1,numodeL)),numcolsL+1);
            diffmeshatswitch(1)*DRal(numbccoord,odecoordL).*psivalL(1,ll,numcolsL+1);
            J(gen_c_start:gen_c)=ans(:);
            ROW(gen_c_start:gen_c)=repmat(ROWcounter+numbccoord.',numodeL,1);
            repmat(COLcounter+odecoordL,numbc,1);
            COL(gen_c_start:gen_c)=ans(:);
            COLcounter=COLcounter+numodeL;
            if ~isempty(aecoordL)
                gen_c_start=gen_c+1;
                gen_c=gen_c+numbc*numaeL;
                %DRal(numbccoord,aecoordL).*psival0L(ones(numbc,1),ll(1,ones(1,numaeL)),numcolsL+1);
                DRal(numbccoord,aecoordL).*psival0L(1,ll,numcolsL+1);
                J(gen_c_start:gen_c)=ans(:);
                ROW(gen_c_start:gen_c)=repmat(ROWcounter+numbccoord.',numaeL,1);
                repmat(COLcounter+(1:numaeL),numbc,1);
                COL(gen_c_start:gen_c)=ans(:);
                COLcounter=COLcounter+numaeL;
            end
        end
    end
    gen_c_start=gen_c+1;
    gen_c=gen_c+numbc*numode;
    DRar(numbccoord,odecoord);
    J(gen_c_start:gen_c)=ans(:);
    repmat(ROWcounter+numbccoord.',1,numode);
    ROW(gen_c_start:gen_c)=ans(:);
    repmat(COLcounter+odecoord,numbc,1);
    COL(gen_c_start:gen_c)=ans(:);
    %COLcounter=COLcounter+(rightarcindex(arc)-leftarcindex(arc)-1)*(numcols*numeq+numode);
    COLcounter=COLcounter+(rightarcindex(arc)-leftarcindex(arc)-arc)*(numcols*numeq+numode);
    COLcounter_aL=COLcounter;

    gen_c_start=gen_c+1;
    gen_c=gen_c+numbc*numode;
    DRbl(numbccoord,odecoord);
    J(gen_c_start:gen_c)=ans(:);
    repmat(ROWcounter+numbccoord.',1,numode);
    ROW(gen_c_start:gen_c)=ans(:);
    repmat(COLcounter+odecoord,numbc,1);
    COL(gen_c_start:gen_c)=ans(:);
    COLcounter=COLcounter+numode;
    for ll=numcolscoord
        gen_c_start=gen_c+1;
        gen_c=gen_c+numbc*numode;
        %diffmeshatswitch(3)*DRbl(numbccoord,odecoord).*psival(ones(numbc,1),ll(1,ones(1,numode)),numcols+1);
        diffmeshatswitch(3)*DRbl(numbccoord,odecoord)*psival(1,ll,numcols+1);
        J(gen_c_start:gen_c)=ans(:);
        repmat(ROWcounter+numbccoord.',1,numode);
        ROW(gen_c_start:gen_c)=ans(:);
        repmat(COLcounter+odecoord,numbc,1);
        COL(gen_c_start:gen_c)=ans(:);
        COLcounter=COLcounter+numode;
        if ~isempty(aecoord)
            gen_c_start=gen_c+1;
            gen_c=gen_c+numbc*numae;
            %DRbl(numbccoord,aecoord).*psival0(ones(numbc,1),ll(1,ones(1,numae)),numcols+1);
            DRbl(numbccoord,aecoord)*psival0(1,ll,numcols+1);
            J(gen_c_start:gen_c)=ans(:);
            ROW(gen_c_start:gen_c)=repmat(ROWcounter+numbccoord.',numae,1);
            repmat(COLcounter+(1:numae),numbc,1);
            COL(gen_c_start:gen_c)=ans(:);
            COLcounter=COLcounter+numae;
        end
    end
    if ~isempty(DRbr)
        gen_c_start=gen_c+1;
        gen_c=gen_c+numbc*numodeR;
        DRbr(numbccoord,odecoordR);
        J(gen_c_start:gen_c)=ans(:);
        ROW(gen_c_start:gen_c)=repmat(ROWcounter+numbccoord.',numodeR,1);
        repmat(COLcounter+odecoordR,numbc,1);
        COL(gen_c_start:gen_c)=ans(:);
        COLcounter=COLcounter+numodeR;
        for ll=numcolscoordR
            gen_c_start=gen_c+1;
            gen_c=gen_c+numbc*numodeR;
            %diffmeshatswitch(4)*DRbr(numbccoord,odecoordR).*psivalR(ones(numbc,1),ll(1,ones(1,numodeR)),numcolsR+2);
            diffmeshatswitch(4)*DRbr(numbccoord,odecoordR).*psivalR(1,ll,numcolsR+2);
            J(gen_c_start:gen_c)=ans(:);
            ROW(gen_c_start:gen_c)=repmat(ROWcounter+numbccoord.',numodeR,1);
            repmat(COLcounter+odecoordR,numbc,1);
            COL(gen_c_start:gen_c)=ans(:);
            COLcounter=COLcounter+numodeR;
            if ~isempty(aecoordR)
                gen_c_start=gen_c+1;
                gen_c=gen_c+numbc*numaeR;
                %DRbr(numbccoord,aecoordR).*psival0R(ones(numbc,1),ll(1,ones(1,numaeR)),numcolsR+2);
                DRbr(numbccoord,aecoordR).*psival0R(1,ll,numcolsR+2);
                J(gen_c_start:gen_c)=ans(:);
                ROW(gen_c_start:gen_c)=repmat(ROWcounter+numbccoord.',numaeR,1);
                repmat(COLcounter+(1:numaeR),numbc,1);
                COL(gen_c_start:gen_c)=ans(:);
                COLcounter=COLcounter+numaeR;
            end
        end
    end
    gen_c_start=gen_c+1;
    gen_c=gen_c+numparameter*numbc;
    DRp(numbccoord,1:numparameter);
    J(gen_c_start:gen_c)=ans(:);
    ROW(gen_c_start:gen_c)=repmat(ROWcounter+numbccoord.',numparameter,1);
    repmat(OCMATCONT.HE.coeffcoordmp+(1:numparameter),numbc,1);
    COL(gen_c_start:gen_c)=ans(:);
    COLcounter=COLcounter_aL;
    ROWcounter=ROWcounter+numbc;
    
    
    %----------Start Jacobian entries at collocations and mesh points
    arcindex=OCMATCONT.HE.arcindex(arc);
    psival=domainddata(arcindex).psival;
    numcols=domainddata(arcindex).numcols;
    numcolscoord=domainddata(arcindex).numcolscoord;
    numeq=domainddata(arcindex).numeq;
    numae=domainddata(arcindex).numae;
    numode=domainddata(arcindex).numode;
    odecoord=domainddata(arcindex).odecoord; % coordinate of (first order) odes
    aecoord=domainddata(arcindex).aecoord; % coordinate of (first order) odes
    eqcoord=domainddata(arcindex).eqcoord;
    eqderivative=eye(numode);
    if numae
        eqderivative=blkdiag(eye(numode),zeros(numae));
    end
    fz=zeros(numeq,2,numcols);
    for step=leftarcindex(arc):rightarcindex(arc)-1
        countermeshpts=step-leftarcindex(arc)+1;
        for kk=numcolscoord
            diffmesh_arc=diffmesh(leftarcindex(arc):rightarcindex(arc)-1);
            fz(odecoord,1,numcolscoord)=evaluatePcols(countermeshpts,diffmesh_arc,y(odecoord,1,idx),z(odecoord,numcolscoord,idx),odecoord,psival,numcols);
            fz(odecoord,2,numcolscoord)=z(odecoord,numcolscoord,idx(countermeshpts));
            fz(aecoord,1,numcolscoord)=z(aecoord,numcolscoord,idx(countermeshpts));
            %ii indicates, which component of g_ik will be differentiated with respect
            %to y und z
            [Dg Dpg]=odejac(collocationpts(countermeshpts,kk),fz(:,1,kk),modelpar,freepar,contval,arc);
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
                    J(gen_c_start:gen_c)=Dpg(ii,1:numparameter);
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
