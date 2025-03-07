function Jac=calcresjac_dae(tcol,x,y,z1,freepar,modelpar,dae,bc,ic,daejac,bcjac,bcic,tangent)
global OCMATCONT

if ~isempty(tangent)
    Jac=spalloc(OCMATCONT.MESHDATA.continuationindex-1,OCMATCONT.MESHDATA.continuationindex,OCMATCONT.MESHDATA.JacobianNonzeroEntriesNumber);  % sparse jacobian matrix
    last_cols=zeros(OCMATCONT.MESHDATA.continuationindex-1,OCMATCONT.freeparameternum);   % accumulator for parameter derivatives
    freeparametercoord=1:OCMATCONT.freeparameternum;
else
    Jac=spalloc(OCMATCONT.MESHDATA.continuationindex-1,OCMATCONT.MESHDATA.continuationindex-1,OCMATCONT.MESHDATA.JacobianNonzeroEntriesNumber);  % sparse jacobian matrix
    last_cols=zeros(OCMATCONT.MESHDATA.continuationindex-1,OCMATCONT.freeparameternum-1);   % accumulator for parameter derivatives
    freeparametercoord=1:OCMATCONT.freeparameternum-1;
end
% BC points
xa=x(:,OCMATCONT.MESHDATA.arcpositioncollocationmesh(1,:));
xb=x(:,OCMATCONT.MESHDATA.arcpositioncollocationmesh(2,:));


if isempty(bcjac) % use numerical approximation
    daebcjacfunch=@bcjacnum;
    Fargbc={xa,xb,freepar,modelpar,bc};
else % use analytical Jacobian
    daebcjacfunch=bcjac;
    Fargbc={xa,xb,freepar,modelpar};
end

if isempty(daejac) % use numerical approximation
    daejacfunch=@daejacnum;
    Fargjac={freepar,modelpar,dae};
else % use analytical Jacobian
    daejacfunch=daejac;
    Fargjac={freepar,modelpar};
end

% calculate derivatives of boundary conditions
[dGdxa,dGdxb,dGdpar] = daebcjacfunch(Fargbc{:});
if ~isempty(tangent)
    last_cols(OCMATCONT.bcidx,:) = dGdpar;
else
    last_cols(OCMATCONT.bcidx,:) = dGdpar(:,freeparametercoord);
end
t=tcol(OCMATCONT.MESHDATA.tmeshidx);

for arc = 1:OCMATCONT.arcnumber
    h=diff(t(OCMATCONT.MESHDATA.arcposition(1,arc):OCMATCONT.MESHDATA.arcposition(2,arc)));
    % generate BC block
    Jac(OCMATCONT.bcidx,OCMATCONT.MESHDATA.ycoefficientidx(:,1))=dGdxa(:,OCMATCONT.firstordercoordinate);
    Jac(OCMATCONT.bcidx,OCMATCONT.MESHDATA.BCbBlockidx)=[dGdxb(:,OCMATCONT.firstordercoordinate) myrepmat(dGdxb,1,OCMATCONT.CollocationNumber,OCMATCONT.bcidx(end),OCMATCONT.componentnumber)].*OCMATCONT.MESHDATA.BCbBlock;

    rows=OCMATCONT.bcidx(end)+1;
    cols=1;
    for ii=1:OCMATCONT.MESHDATA.meshNumber(arc)-1
        % generate J block
        collidx=((ii-1)*(OCMATCONT.CollocationNumber+1)+1):(ii*(OCMATCONT.CollocationNumber+1)+1);
        tcoll=tcol(collidx);
        xcoll=x(:,collidx);
        absoluteidx=rel2absidx(OCMATCONT.BasicJBlockSize,[rows,cols],OCMATCONT.MESHDATA.continuationindex-1);
        JTl=zeros(OCMATCONT.CollocationNumber*OCMATCONT.componentnumber,OCMATCONT.componentnumber); % Jacobians at collocation points Til
        for jj=1:OCMATCONT.CollocationNumber
            [tmp,Jpar]=daejacfunch(tcoll(jj+1),xcoll(:,jj+1),Fargjac{:});
            tmp(OCMATCONT.firstordercoordinate,:)=-tmp(OCMATCONT.firstordercoordinate,:);
            JTl((jj-1)*OCMATCONT.componentnumber+1:jj*OCMATCONT.componentnumber,:)=tmp;
            Jpar(OCMATCONT.firstordercoordinate,:)=-Jpar(OCMATCONT.firstordercoordinate,:);
            last_cols(rows:rows+OCMATCONT.componentnumber-1,:)=Jpar(:,freeparametercoord);
            rows=rows+OCMATCONT.componentnumber;
        end
        BasicJBlock=OCMATCONT.BasicJBlock;
        BasicJBlock(:,OCMATCONT.JBlockidx)=h(ii)*OCMATCONT.JBlock;
        JBlock0=[JTl(:,OCMATCONT.firstordercoordinate,:) myrepmat(JTl,1,OCMATCONT.CollocationNumber,OCMATCONT.CollocationNumber*OCMATCONT.componentnumber,OCMATCONT.componentnumber)];
        JBlock0=BasicJBlock.*JBlock0;
        JBlock0(OCMATCONT.BasicJBlockAddOneidx)=JBlock0(OCMATCONT.BasicJBlockAddOneidx)+1;
        Jac(absoluteidx(OCMATCONT.BasicJBlocknonzeroidx))=JBlock0(OCMATCONT.BasicJBlocknonzeroidx);

        if ii<OCMATCONT.MESHDATA.meshNumber(arc)-1
            %generate C block
            BasicCBlock=OCMATCONT.BasicCBlock;
            CBlock=OCMATCONT.CBlock*h(ii);
            BasicCBlock(:,OCMATCONT.JBlockidx)=CBlock;
            absoluteidx=rel2absidx(OCMATCONT.BasicCBlockSize,[rows,cols],OCMATCONT.MESHDATA.continuationindex-1);
            Jac(absoluteidx(OCMATCONT.BasicCBlocknonzeroidx))=BasicCBlock(OCMATCONT.BasicCBlocknonzeroidx);
        end

        rows=rows+OCMATCONT.firstordercomponentnumber;
        cols=cols+OCMATCONT.firstordercomponentnumber+OCMATCONT.CollocationNumber*OCMATCONT.componentnumber;
    end
end
% add parameter derivatives
if ~isempty(tangent)
    Jac(:,cols:cols+OCMATCONT.freeparameternum-1)=last_cols;
    Jac(OCMATCONT.MESHDATA.continuationindex,:)=tangent(:).';
else
    Jac(:,cols:cols+OCMATCONT.freeparameternum-2)=last_cols(:,freeparametercoord);
end

function [J,Jpar]=daejacnum(t,x,freepar,modelpar,dae)
global OCMATCONT

numJacOpt.diffvar=2;
numJacOpt.vectvars=[];
J=numjaccsd(dae,{t,x,freepar,modelpar},OCMATCONT.componentnumber,numJacOpt);

numJacOpt.diffvar=3;
Jpar=numjaccsd(dae,{t,x,freepar,modelpar},OCMATCONT.componentnumber,numJacOpt);


function [dGdxa,dGdxb,dGdpar]=bcjacnum(xa,xb,freepar,modelpar,bc)
global OCMATCONT

numJacOpt.diffvar=1;
numJacOpt.vectvars=[];
dGdxa=numjaccsd(bc,{xa,xb,freepar,modelpar},OCMATCONT.firstordercomponentnumber,numJacOpt);

numJacOpt.diffvar=2;
dGdxb=numjaccsd(bc,{xa,xb,freepar,modelpar},OCMATCONT.firstordercomponentnumber,numJacOpt);

numJacOpt.diffvar=3;
dGdpar=numjaccsd(bc,{xa,xb,freepar,modelpar},OCMATCONT.firstordercomponentnumber,numJacOpt);

