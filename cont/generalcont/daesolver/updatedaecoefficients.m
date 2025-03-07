function updatedaecoefficients(tmesh,substructname)
global OCMATCONT
% returns the coefficients for MESHDATA1_2 correspoing to coeff for
% MESHDATA
ctr=0;
ctr1=0;
ctr2=0;
ctr3=0;
addidx=0;
tst=0;
% positional indices for the equations
% bcnum ... number of boundary conditions
% firstorderFinterioridx
% zeroorderFinterioridx
% transitionFidx
% ctrfoFi
% ctrzoFi
% ctrtF

ctrfoFi=0;
ctrzoFi=0;
ctrtF=0;

[dum,tmeshidx]=makecollocationmesh(tmesh);
OCMATCONT.(substructname).tmeshidx=tmeshidx;
idx=find(diff(tmesh)==0);

N=length(tmesh);
arcposition=[1 idx+1;idx N];
OCMATCONT.(substructname).arcposition=arcposition;

N=arcposition(2,:);
OCMATCONT.(substructname).meshNumber=N;
idx=cumsum((N-1)*(OCMATCONT.CollocationNumber+1)+1);
OCMATCONT.(substructname).arcpositioncollocationmesh=[1 idx(1:OCMATCONT.arcnumber-1)+1;idx];


coeff=zeros(sum(N-1)*((OCMATCONT.CollocationNumber+1)*OCMATCONT.firstordercomponentnumber+OCMATCONT.CollocationNumber*OCMATCONT.zeroordercomponentnumber)+OCMATCONT.freeparameternum,1);

tmp=arcposition(2,:)-(1:OCMATCONT.arcnumber);
OCMATCONT.(substructname).diffarcposition=[1 tmp(1:end-1)+1;tmp];

firstorderidx=zeros(OCMATCONT.CollocationNumber+1,OCMATCONT.firstordercomponentnumber*OCMATCONT.(substructname).diffarcposition(2,OCMATCONT.arcnumber));
zeroorderidx=zeros(OCMATCONT.CollocationNumber,OCMATCONT.zeroordercomponentnumber*OCMATCONT.(substructname).diffarcposition(2,OCMATCONT.arcnumber));

addFpos=OCMATCONT.bcidx(end);
OCMATCONT.(substructname).firstorderFinterioridx=zeros(OCMATCONT.firstordercomponentnumber*OCMATCONT.CollocationNumber*sum(N-1),1);
OCMATCONT.(substructname).zeroorderFinterioridx=zeros(OCMATCONT.zeroordercomponentnumber*OCMATCONT.CollocationNumber*sum(N-1),1);
OCMATCONT.(substructname).transitionFidx=zeros(OCMATCONT.firstordercomponentnumber*sum(N-2),1);

arcpart4firstorderFinterioridx=OCMATCONT.firstordercomponentnumber*OCMATCONT.CollocationNumber*cumsum(N-1);
OCMATCONT.(substructname).arcpart4firstorderFinterioridx=[1 arcpart4firstorderFinterioridx(1:end-1)+1;arcpart4firstorderFinterioridx];

arcpart4zeroorderFinterioridx=OCMATCONT.zeroordercomponentnumber*OCMATCONT.CollocationNumber*cumsum(N-1);
OCMATCONT.(substructname).arcpart4zeroorderFinterioridx=[1 arcpart4zeroorderFinterioridx(1:end-1)+1;arcpart4zeroorderFinterioridx];

arcpart4transitionFidx=OCMATCONT.firstordercomponentnumber*cumsum(N-2);
OCMATCONT.(substructname).arcpart4transitionFidx=[1 arcpart4transitionFidx(1:end-1)+1;arcpart4transitionFidx];

for arc=1:OCMATCONT.arcnumber
    % equation index, total number of equations for first order ODEs and
    % AEs for each arc
    eqnum=OCMATCONT.componentnumber*OCMATCONT.CollocationNumber*(N(arc)-1)+OCMATCONT.firstordercomponentnumber*(N(arc)-2);
    Fidx=(1:eqnum);
    
    % equations at the mesh points
    transitionFidx=repmat([zeros(OCMATCONT.componentnumber*OCMATCONT.CollocationNumber,1);ones(OCMATCONT.firstordercomponentnumber,1)],N(arc)-1,1);
    transitionFidx=transitionFidx(1:eqnum);
    transitionFidx=find(transitionFidx==1);
    Fidx(transitionFidx)=[];
    Fidx=reshape(addFpos+Fidx,OCMATCONT.componentnumber,[]);
    OCMATCONT.(substructname).transitionFidx(ctrtF+(1:OCMATCONT.firstordercomponentnumber*(N(arc)-2)))=addFpos+transitionFidx;
    OCMATCONT.(substructname).firstorderFinterioridx(ctrfoFi+(1:OCMATCONT.firstordercomponentnumber*OCMATCONT.CollocationNumber*(N(arc)-1)))=Fidx(OCMATCONT.firstordercoordinate,:);
    OCMATCONT.(substructname).zeroorderFinterioridx(ctrzoFi+(1:OCMATCONT.zeroordercomponentnumber*OCMATCONT.CollocationNumber*(N(arc)-1)))=Fidx(OCMATCONT.zeroordercoordinate,:);
    addFpos=addFpos+eqnum;
    ctrfoFi=ctrfoFi+OCMATCONT.firstordercomponentnumber*OCMATCONT.CollocationNumber*(N(arc)-1);
    ctrzoFi=ctrzoFi+OCMATCONT.zeroordercomponentnumber*OCMATCONT.CollocationNumber*(N(arc)-1);
    ctrtF=ctrtF+OCMATCONT.firstordercomponentnumber*(N(arc)-2);
    
    tst=tst+(OCMATCONT.CollocationNumber+1)*(N(arc)-1);

    
    locctr=0;
    locctr1=1;
    coeff0=NaN(OCMATCONT.zeroordercomponentnumber,(OCMATCONT.CollocationNumber+1)*(N(arc)-1));
    for ii=1:N(arc)-1
        coeff0(:,locctr1+(1:(OCMATCONT.CollocationNumber)))=0;
        locctr1=locctr1+OCMATCONT.CollocationNumber+1;
        locctr=locctr+OCMATCONT.CollocationNumber-1;
    end

    locctr=0;
    locctr1=0;
    coeff1=zeros(OCMATCONT.firstordercomponentnumber,(OCMATCONT.CollocationNumber+1)*(N(arc)-1));
    for ii=1:N(arc)-1
        coeff1(:,locctr1+(1:(OCMATCONT.CollocationNumber+1)))=1;
        locctr=locctr+OCMATCONT.CollocationNumber;
        locctr1=locctr1+OCMATCONT.CollocationNumber+1;
    end
    coeffarc(OCMATCONT.firstordercoordinate,:)=coeff1;
    coeffarc(OCMATCONT.zeroordercoordinate,:)=coeff0;
    coeffarc=coeffarc(:);
    idx=isnan(coeffarc);
    coeffarc(idx)=[];

    coeff(ctr+(1:(N(arc)-1)*((OCMATCONT.CollocationNumber+1)*OCMATCONT.firstordercomponentnumber+OCMATCONT.CollocationNumber*OCMATCONT.zeroordercomponentnumber)))=coeffarc;
    
    ctr=ctr+(N(arc)-1)*((OCMATCONT.CollocationNumber+1)*OCMATCONT.firstordercomponentnumber+OCMATCONT.CollocationNumber*OCMATCONT.zeroordercomponentnumber);

    ctr1=ctr1+(N(arc)-1)*(OCMATCONT.CollocationNumber+1)+1;


    coeff0=ones(OCMATCONT.zeroordercomponentnumber,(OCMATCONT.CollocationNumber+1)*(N(arc)-1));
    coeff0(:,1:OCMATCONT.CollocationNumber+1:end)=0;
    coeff1=ones(OCMATCONT.firstordercomponentnumber,(OCMATCONT.CollocationNumber+1)*(N(arc)-1));

    coeffidx(OCMATCONT.firstordercoordinate,:)=coeff1;
    coeffidx(OCMATCONT.zeroordercoordinate,:)=coeff0;
    coeffidx=cumsum(coeffidx(:));
    coeffidx=reshape(coeffidx,OCMATCONT.componentnumber,[]);
    coeffidx(idx)=NaN;
    coeffidx=coeffidx+addidx;
    % index to retrieve first and zero order components from the coefficients
    firstorderidx(:,ctr2+(1:OCMATCONT.firstordercomponentnumber*(N(arc)-1)))=reshape(permute(reshape(coeffidx(OCMATCONT.firstordercoordinate,:),OCMATCONT.firstordercomponentnumber,OCMATCONT.CollocationNumber+1,N(arc)-1),[2 1 3]),OCMATCONT.CollocationNumber+1,[]);
    tmpzeroorderidx=reshape(permute(reshape(coeffidx(OCMATCONT.zeroordercoordinate,:),OCMATCONT.zeroordercomponentnumber,OCMATCONT.CollocationNumber+1,N(arc)-1),[2 1 3]),OCMATCONT.CollocationNumber+1,[]);
    tmpzeroorderidx(isnan(tmpzeroorderidx))=[];
    tmpzeroorderidx=reshape(tmpzeroorderidx,OCMATCONT.CollocationNumber,[]);
    zeroorderidx(:,ctr3+(1:OCMATCONT.zeroordercomponentnumber*(N(arc)-1)))=tmpzeroorderidx;
    addidx=addidx+max([firstorderidx(:);zeroorderidx(:)]);
    ctr2=ctr2+OCMATCONT.firstordercomponentnumber*(N(arc)-1);
    ctr3=ctr3+OCMATCONT.zeroordercomponentnumber*(N(arc)-1);
end
% 'firstorderidx' and 'zeororderidx' retrieve the coefficients of the order
% one and zero from the total coeffcient vector in the form (n ... number
% of first order components, m ... number of zero order components, N mesh
% size:
% firstorderidx: (tau_i T_{i,1} ... T_{i,m}) x n*(N-1) 
% zeroorderidx: (T_{i,1} ... T_{i,m}) x n*(N-1) 
OCMATCONT.(substructname).firstorderidx=firstorderidx;
OCMATCONT.(substructname).zeroorderidx=zeroorderidx;

% the arc partition for the 'firstorderidx' and 'zeororderidx'
arcpart4firstorderidx=cumsum(OCMATCONT.firstordercomponentnumber*(N(arc)-1));
OCMATCONT.(substructname).arcpart4firstorderidx=[1 arcpart4firstorderidx(1:end-1)+1;arcpart4firstorderidx];
arcpart4zeroorderidx=cumsum(OCMATCONT.zeroordercomponentnumber*(N(arc)-1));
OCMATCONT.(substructname).arcpart4zeroorderidx=[1 arcpart4zeroorderidx(1:end-1)+1;arcpart4zeroorderidx];

% the y and first order z coefficients are directly used for the system of
% equations. The according coefficients allow to retrieve these
% coefficients immediately from the total coefficient vector in the matrix
% form (k ... number collocation points): 
% ycoefficientidx : n x N-1
% firstorderzcoefficientidx : n x k*(N-1)
OCMATCONT.(substructname).ycoefficientidx=reshape(OCMATCONT.(substructname).firstorderidx(1,:),OCMATCONT.firstordercomponentnumber,[]); % index to retrieve the y coefficients
OCMATCONT.(substructname).firstorderzcoefficientidx=reshape(permute(reshape(OCMATCONT.(substructname).firstorderidx(2:OCMATCONT.CollocationNumber+1,:),OCMATCONT.CollocationNumber,OCMATCONT.firstordercomponentnumber,sum(N-1)),[2 1 3]),OCMATCONT.firstordercomponentnumber,[]); % index to retrieve the z coefficients for first order components

% the arc partition for the 'ycoefficientidx' and 'firstorderzcoefficientidx'
arcpart4ycoefficientidx=cumsum(N(arc)-1);
OCMATCONT.(substructname).arcpart4ycoefficientidx=[1 arcpart4ycoefficientidx(1:end-1)+1;arcpart4ycoefficientidx];

arcpart4firstorderzcoefficientidx=cumsum(OCMATCONT.CollocationNumber*(N-1));
OCMATCONT.(substructname).arcpart4firstorderzcoefficientidx=[1 arcpart4firstorderzcoefficientidx(1:end-1)+1;arcpart4firstorderzcoefficientidx];


OCMATCONT.(substructname).freeparameterindex=ctr+(1:OCMATCONT.freeparameternum); % position of the free parameters in the coefficient vector
OCMATCONT.(substructname).continuationindex=ctr+OCMATCONT.freeparameternum;

% JBlock ... block of psi entries for first order components derived from
% the derivative with respect to c_ikj
% BasicJBlock ... ones derived from the derivative with respect to y_ik,
% JBlockidx ... is position index of JBlock in BasicJBlock, i.e.
% BasicJBlock(JBlockidx)=JBlock 
% 
% BasicJBlocknonzeroidx ... indices of nonzero elements of BasicJBlock
% BasicJBlockAddOneidx ... indices where a one is added to the derivative
% of F_r(Til) with respect to c_iks, since r=k and l=s.
% OCMATCONT.(substructname).JBlock=zeros(OCMATCONT.componentnumber*OCMATCONT.CollocationNumber,OCMATCONT.firstordercomponentnumber*OCMATCONT.CollocationNumber);
% OCMATCONT.(substructname).JBlockidx=OCMATCONT.firstordercomponentnumber+OCMATCONT.firstordercoordinate(:,ones(1,OCMATCONT.CollocationNumber))+OCMATCONT.componentnumber*repmat(0:OCMATCONT.CollocationNumber-1,OCMATCONT.firstordercomponentnumber,1);
% OCMATCONT.(substructname).JBlockidx=OCMATCONT.(substructname).JBlockidx(:).';
% OCMATCONT.(substructname).BasicJBlock=zeros(OCMATCONT.componentnumber*OCMATCONT.CollocationNumber,OCMATCONT.CollocationNumber*OCMATCONT.componentnumber+OCMATCONT.firstordercomponentnumber);
% OCMATCONT.(substructname).BasicJBlock(:,1:OCMATCONT.firstordercomponentnumber)=1;
% 
% OCMATCONT.(substructname).BasicJBlocknonzeroidx=zeros(1,OCMATCONT.zeroordercomponentnumber*OCMATCONT.componentnumber*(OCMATCONT.CollocationNumber-1));
% OCMATCONT.(substructname).BasicJBlockAddOneidx=zeros(1,OCMATCONT.firstordercomponentnumber*OCMATCONT.CollocationNumber);
% OCMATCONT.(substructname).BasicJBlockSize=[OCMATCONT.CollocationNumber*OCMATCONT.componentnumber,OCMATCONT.CollocationNumber*OCMATCONT.componentnumber+OCMATCONT.firstordercomponentnumber];
% 
% OCMATCONT.(substructname).CBlock=zeros(OCMATCONT.firstordercomponentnumber,OCMATCONT.firstordercomponentnumber*OCMATCONT.CollocationNumber);
% OCMATCONT.(substructname).BasicCBlock=zeros(OCMATCONT.firstordercomponentnumber,OCMATCONT.CollocationNumber*OCMATCONT.componentnumber+2*OCMATCONT.firstordercomponentnumber);
% 
% OCMATCONT.(substructname).BasicCBlock(:,1:OCMATCONT.firstordercomponentnumber)=eye(OCMATCONT.firstordercomponentnumber);
% OCMATCONT.(substructname).BasicCBlock(:,OCMATCONT.CollocationNumber*OCMATCONT.componentnumber+OCMATCONT.firstordercomponentnumber+(1:OCMATCONT.firstordercomponentnumber))=-eye(OCMATCONT.firstordercomponentnumber);
% OCMATCONT.(substructname).BasicCBlockSize=[OCMATCONT.firstordercomponentnumber,OCMATCONT.CollocationNumber*OCMATCONT.componentnumber+2*OCMATCONT.firstordercomponentnumber];

OCMATCONT.(substructname).BCbBlock=ones(OCMATCONT.firstordercomponentnumber+OCMATCONT.freeparameternum-1,OCMATCONT.firstordercomponentnumber+OCMATCONT.componentnumber*OCMATCONT.CollocationNumber);
OCMATCONT.(substructname).BCbBlockidx=zeros(1,OCMATCONT.firstordercomponentnumber+OCMATCONT.componentnumber*OCMATCONT.CollocationNumber);

%tic
%psival=OCMATCONT.FirstOrderCollocationValues(2:OCMATCONT.CollocationNumber+1,2:OCMATCONT.CollocationNumber+1);
%psival=psival.';
%psival=psival(:).';
%psival=psival(ones(OCMATCONT.componentnumber*OCMATCONT.firstordercomponentnumber,1),:);

ctr=0;
% rows=1:OCMATCONT.componentnumber;
% cols=1:OCMATCONT.firstordercomponentnumber;
% tmp=zeros(OCMATCONT.componentnumber,OCMATCONT.firstordercomponentnumber);
% tmp1=ones(OCMATCONT.componentnumber,OCMATCONT.zeroordercomponentnumber);
% tmp2=cumsum([0 (OCMATCONT.firstordercoordinate(2:OCMATCONT.firstordercomponentnumber)-OCMATCONT.firstordercoordinate(1:OCMATCONT.firstordercomponentnumber-1)).']*OCMATCONT.componentnumber*OCMATCONT.CollocationNumber);

tmpb1=zeros(OCMATCONT.firstordercomponentnumber);
tmpb0=zeros(OCMATCONT.firstordercomponentnumber,OCMATCONT.zeroordercomponentnumber);

h=diff(tmesh(arcposition(1,1):arcposition(2,1)));

for ii=0:OCMATCONT.CollocationNumber-1
    for jj=0:OCMATCONT.CollocationNumber-1
        ctr=ctr+1;
        %tmp(:)=psival(:,ctr);
        %OCMATCONT.(substructname).JBlock(ii*OCMATCONT.componentnumber+rows,jj*OCMATCONT.firstordercomponentnumber+cols)=tmp;
        %if ii==jj
        %    OCMATCONT.(substructname).BasicJBlock(ii*OCMATCONT.componentnumber+rows,OCMATCONT.firstordercomponentnumber+jj*OCMATCONT.componentnumber+OCMATCONT.zeroordercoordinate)=tmp1;
        %end
    end
    %OCMATCONT.(substructname).BasicJBlockAddOneidx(ii*OCMATCONT.firstordercomponentnumber+1:(ii+1)*OCMATCONT.firstordercomponentnumber)=(OCMATCONT.firstordercomponentnumber+ii*OCMATCONT.componentnumber)*OCMATCONT.componentnumber*OCMATCONT.CollocationNumber+ii*OCMATCONT.componentnumber+OCMATCONT.firstordercoordinate(:).'+tmp2;
    
    %OCMATCONT.(substructname).CBlock(:,ii*OCMATCONT.firstordercomponentnumber+cols)=eye(OCMATCONT.firstordercomponentnumber)*OCMATCONT.FirstOrderCollocationValues(ii+2,OCMATCONT.CollocationNumber+2);
    
    tmpb1(:)=h(end-1)*OCMATCONT.FirstOrderCollocationValues(ii+2,OCMATCONT.CollocationNumber+2);
    tmpb0(:)=OCMATCONT.ZeroOrderCollocationValues(ii+1,OCMATCONT.CollocationNumber+2);
    OCMATCONT.(substructname).BCbBlock(:,OCMATCONT.firstordercomponentnumber+ii*OCMATCONT.componentnumber+OCMATCONT.firstordercoordinate)=tmpb1;
    OCMATCONT.(substructname).BCbBlock(:,OCMATCONT.firstordercomponentnumber+ii*OCMATCONT.componentnumber+OCMATCONT.zeroordercoordinate)=tmpb0;
end

% BasicJBlock=OCMATCONT.(substructname).BasicJBlock;
% BasicCBlock=OCMATCONT.(substructname).BasicCBlock;
% BasicJBlock(:,OCMATCONT.(substructname).JBlockidx)=OCMATCONT.(substructname).JBlock;
% BasicJBlock(OCMATCONT.(substructname).BasicJBlockAddOneidx)=1;
% BasicCBlock(:,OCMATCONT.(substructname).JBlockidx)=OCMATCONT.(substructname).CBlock;
% OCMATCONT.(substructname).BasicJBlocknonzeroidx=find(BasicJBlock).';
% %OCMATCONT.BasicJBlockzeroidx=find(BasicJBlock==0).';
% OCMATCONT.(substructname).BasicCBlocknonzeroidx=find(BasicCBlock).';
%OCMATCONT.BasicCBlockzeroidx=find(BasicCBlock==0).';

OCMATCONT.(substructname).BCbBlockidx=(N-2)*(OCMATCONT.componentnumber*OCMATCONT.CollocationNumber+OCMATCONT.firstordercomponentnumber)+ ...
    (1:OCMATCONT.firstordercomponentnumber+OCMATCONT.componentnumber*OCMATCONT.CollocationNumber);

OCMATCONT.(substructname).JacobianNonzeroEntriesNumber=(N-1)*OCMATCONT.componentnumber*OCMATCONT.CollocationNumber*(OCMATCONT.firstordercomponentnumber*OCMATCONT.CollocationNumber+1)+ ...
    (N-2)*OCMATCONT.firstordercomponentnumber*(OCMATCONT.CollocationNumber+2)+ ...
    (OCMATCONT.firstordercomponentnumber+OCMATCONT.freeparameternum-1)*(OCMATCONT.componentnumber*OCMATCONT.CollocationNumber+ ...
    2*OCMATCONT.firstordercomponentnumber)+((N-1)*OCMATCONT.componentnumber*OCMATCONT.CollocationNumber+OCMATCONT.firstordercomponentnumber+ ...
    OCMATCONT.freeparameternum-1)*OCMATCONT.freeparameternum;

if strcmp(substructname,'MESHDATA1_2')
    for ii=1:OCMATCONT.arcnumber
        keepidx=2*(0:((OCMATCONT.MESHDATA1_2.meshNumber(ii)-2)/2))*(OCMATCONT.NumberPoints2_1-1);
        keepidx=[repmat(keepidx(:),1,numel(OCMATCONT.IndexLeft2_1))+repmat(OCMATCONT.IndexLeft2_1,length(keepidx),1) repmat(keepidx(:),1,numel(OCMATCONT.IndexRight2_1)-1)+OCMATCONT.NumberPoints2_1-1+repmat(OCMATCONT.IndexRight2_1(1:end-1),length(keepidx),1)];
        keepidx=keepidx.';
        OCMATCONT.MESHDATA1_2.Matrix(ii).KeepIdx2_1=[keepidx(:).' (OCMATCONT.MESHDATA1_2.meshNumber(ii)-1)*(OCMATCONT.NumberPoints2_1-1)+1];
    end
end