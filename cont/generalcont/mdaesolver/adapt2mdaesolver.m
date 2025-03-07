function sol=adapt2mdaesolver(sol,tnew,order,collMethod,collocnum,interpMethod)
% 
% ADAPT2MDAESOLVER uses the solution specified in 'sol' to initialize the
% Runge-Kutta coeffecients at the time mesh 'tnew'. If tnew=[] the
% coefficients are computed for the time mesh specified in 'sol.x'.
%
% Cells characterize different stages, i.e. sol.x{1},...,sol.x{ns} and
% sol.y{1},...,sol.y{ns}. Each stage can consist of different arcs, these
% are charactized by the time mesh having the form, t^i=sol.x{i} and being
% normalized to [0,...,1],[1,...,2],...,[k-1,...,k], where k is the number
% of arcs. If the time mesh is not normalized it will be normalized during
% the initilization process.
%
% Every variable that can change at different changes and in general cannot
% be represented by a matrix/vector is formulated as a cell.
global OCMATCONT

if isempty(sol) || ~isfield(sol,'x') || ~isfield(sol,'y')
    sol=[];
    return
end
if ~isecell(sol.x)
    t{1}=sol.x;
    X{1}=sol.y;
    stagenum=1;
else
    t=sol.x;
    X=sol.y;
    stagenum=length(t);
end

if isempty(t) || isempty(X)
    sol=[];
    return
end

for ii=1:stagenum
    if length(t{ii})~=size(X{ii},2)
        sol=[];
        return
    end
end
if ~iscell(order)
    ordertmp=order;
    clear order
    order{1}=ordertmp;
    order=order(ones(1,stagenum));
end
if stagenum~=length(order)
    sol=[];
    return
end

if isempty(tnew)
    tnew=t;
end

if isfield(sol,'parameters')
    freepar=sol.parameters;
else
    freepar=[];
    sol.parameters=[];
end

if isempty(freepar)
    ocmatmsg('No free parameter specified. Set dummy parameter.\n')
    freepar=0;
    sol.parameters=0;
end

if nargin<6
    opt=defaultocoptions;
    interpMethod=getocoptions(opt,'SBVPOC','InterpolationMethod');
end

if nargin<5
    collocnum=getocoptions(opt,'SBVPOC','CollocationNumber');
end

if nargin<4
    collMethod=getocoptions(opt,'SBVPOC','CollocationMethod');
end

if isempty(collMethod)
    opt=defaultocoptions;
    collMethod=getocoptions(opt,'SBVPOC','CollocationMethod');
end

if isempty(collocnum)
    opt=defaultocoptions;
    collocnum=getocoptions(opt,'SBVPOC','CollocationNumber');
end

if isempty(interpMethod)
    opt=defaultocoptions;
    interpMethod=getocoptions(opt,'SBVPOC','InterpolationMethod');
end

OCMATCONT.freeparameternum=length(freepar);

OCMATCONT.order=order;
for ii=1:stagenum
    OCMATCONT.firstordercoordinate{ii}=find(order{ii}==1);
    OCMATCONT.zeroordercoordinate{ii}=find(order{ii}==0);
    OCMATCONT.componentnumber(ii)=length(order{ii});
    OCMATCONT.firstordercomponentnumber(ii)=length(OCMATCONT.firstordercoordinate{ii});
    OCMATCONT.zeroordercomponentnumber(ii)=length(OCMATCONT.zeroordercoordinate{ii});
end
OCMATCONT.CollocationNumber=collocnum;
OCMATCONT.CollocationMethod=collMethod;
OCMATCONT.InterpolationMethod=interpMethod;

init_collocationdata();

sol.data.xcol=init_mesh(tnew);

coeff=init_coefficients(sol);
if isfield(sol.data,'coeff') && length(sol.data.coeff)==length(coeff)
    % use existing coeff
    coeff=sol.data.coeff;
end

init_jacobian(sol);

% calculate MESHDATA for fine time grid
% Add a mesh point in-between each two mesh point
sol1_2.x=sol.x;
sol1_2.y=sol.y;
sol1_2.parameters=sol.parameters;

t1_2=makefinemesh(sol.xnew);
idx=find(diff(t1_2)==0);
sol1_2.arcposition=[1 idx+1;idx length(t1_2)];
sol1_2.xnew=t1_2;

[sol1_2,tmeshidx1_2]=init_mesh(sol1_2);
OCMATCONT.MESHDATA1_2.tmeshidx=tmeshidx1_2;

[coeff1_2,tfine1_2]=init_coefficients(sol1_2);
OCMATCONT.MESHDATA1_2.collocationtmesh=tfine1_2;
OCMATCONT.MESHDATA1_2.tmesh=sol1_2.xnew;

init_jacobian(sol);

sol.x=sol.xnew;
sol.y=coeff2points(xcol,coeff,'grid');
sol=rmfield(sol,'xnew');
sol.data.coeff=coeff;
sol.data.xcol=xcol;


%

for ii=1:OCMATCONT.arcnumber
    keepidx=2*(0:((OCMATCONT.MESHDATA1_2.meshNumber(ii)-2)/2))*(OCMATCONT.NumberPoints2_1-1);
    keepidx=[repmat(keepidx(:),1,numel(OCMATCONT.IndexLeft2_1))+repmat(OCMATCONT.IndexLeft2_1,length(keepidx),1) repmat(keepidx(:),1,numel(OCMATCONT.IndexRight2_1)-1)+OCMATCONT.NumberPoints2_1-1+repmat(OCMATCONT.IndexRight2_1(1:end-1),length(keepidx),1)];
    keepidx=keepidx.';
    OCMATCONT.MESHDATA1_2.Matrix(ii).KeepIdx2_1=[keepidx(:).' (OCMATCONT.MESHDATA1_2.meshNumber(ii)-1)*(OCMATCONT.NumberPoints2_1-1)+1];
end


function init_collocationdata()
global OCMATCONT


rho=getcollocationpoints(OCMATCONT.CollocationMethod,OCMATCONT.CollocationNumber);
psi=zeros(OCMATCONT.CollocationNumber,1+OCMATCONT.CollocationNumber);
for ii=1:OCMATCONT.CollocationNumber
    psi(ii,1:OCMATCONT.CollocationNumber+1)=lagrangepolynomials(ii,rho,1);
end

rho_l=rho;
rho_l(rho>=0.5)=[];
rho_g=rho;
rho_g(rho<0.5)=[];
rho_l=[0 rho_l/0.5];
rho_g=[(rho_g-0.5)/0.5 1];
rho_t=[rho_l rho_g];

psival=zeros(OCMATCONT.CollocationNumber+1,OCMATCONT.CollocationNumber+2);
psival1_2=zeros(OCMATCONT.CollocationNumber+1,2*(OCMATCONT.CollocationNumber+1)+1);
psival2_1=zeros(OCMATCONT.CollocationNumber+1,numel(rho_t));
psival(1,:)=1; % first row corresponds to values of the polynomial derived by the integration of psi.   
psival1_2(1,:)=1; % polynomial values for refined grid
psival2_1(1,:)=1; % polynomial values for coarse grid
for ii=1:OCMATCONT.CollocationNumber
    %evaluation of psi at the collocation points,0 and 1 for order 1
    %components
    psival(ii+1,1+(1:OCMATCONT.CollocationNumber))=polyval(psi(ii,:),rho(1:OCMATCONT.CollocationNumber));
    psival(ii+1,OCMATCONT.CollocationNumber+2)=polyval(psi(ii,:),1);

    psival1_2(ii+1,1+(1:OCMATCONT.CollocationNumber))=polyval(psi(ii,:),0.5*rho(1:OCMATCONT.CollocationNumber));
    psival1_2(ii+1,OCMATCONT.CollocationNumber+2+(1:OCMATCONT.CollocationNumber))=polyval(psi(ii,:),0.5*(1+rho(1:OCMATCONT.CollocationNumber)));
    psival1_2(ii+1,OCMATCONT.CollocationNumber+2)=polyval(psi(ii,:),0.5);
    psival1_2(ii+1,2*(OCMATCONT.CollocationNumber+1)+1)=psival(ii+1,OCMATCONT.CollocationNumber+2);

    psival2_1(ii+1,:)=polyval(psi(ii,:),rho_t);
end

%evaluation of psi at the collocation points,0 and 1 for order 0
%components
psi0=zeros(OCMATCONT.CollocationNumber,OCMATCONT.CollocationNumber);
for ii=1:OCMATCONT.CollocationNumber
    psi0(ii,:)=lagrangepolynomials(ii,rho,0);
end
psi0val=zeros(OCMATCONT.CollocationNumber,OCMATCONT.CollocationNumber+2);
psi0val1_2=zeros(OCMATCONT.CollocationNumber,2*(OCMATCONT.CollocationNumber+1)+1);
psi0val2_1=zeros(OCMATCONT.CollocationNumber,numel(rho_t));
for ii=1:OCMATCONT.CollocationNumber
    %evaluation of psi at the collocation points, 0 and 1 for order 0
    %components
    psi0val(ii,1)=polyval(psi0(ii,:),0);
    psi0val(ii,1+(1:OCMATCONT.CollocationNumber))=polyval(psi0(ii,:),rho(1:OCMATCONT.CollocationNumber));
    psi0val(ii,OCMATCONT.CollocationNumber+2)=polyval(psi0(ii,:),1);

    psi0val1_2(ii,1+(1:OCMATCONT.CollocationNumber))=polyval(psi0(ii,:),0.5*rho(1:OCMATCONT.CollocationNumber));
    psi0val1_2(ii,OCMATCONT.CollocationNumber+2+(1:OCMATCONT.CollocationNumber))=polyval(psi0(ii,:),0.5*(1+rho(1:OCMATCONT.CollocationNumber)));
    psi0val1_2(ii,OCMATCONT.CollocationNumber+2)=polyval(psi0(ii,:),0.5);
    psi0val1_2(ii,1)=psi0val(ii,1);
    psi0val1_2(ii,2*(OCMATCONT.CollocationNumber+1)+1)=psi0val(ii,OCMATCONT.CollocationNumber+2);

    psi0val2_1(ii,:)=polyval(psi0(ii,:),rho_t);
end

OCMATCONT.CollocationPoint=rho;
OCMATCONT.FirstOrderCollocationPolynomial=psi;
OCMATCONT.ZeroOrderCollocationPolynomial=psi0;
OCMATCONT.FirstOrderCollocationValues=psival;
OCMATCONT.ZeroOrderCollocationValues=psi0val;
OCMATCONT.FirstOrderCollocationValues1_2=psival1_2;
OCMATCONT.ZeroOrderCollocationValues1_2=psi0val1_2;
OCMATCONT.FirstOrderCollocationValues2_1=psival2_1;
OCMATCONT.ZeroOrderCollocationValues2_1=psi0val2_1;
% OCMATCONT.IndexLeft2_1=1:length(rho_l);
% OCMATCONT.IndexRight2_1=length(rho_l)+(1:length(rho_g));
% OCMATCONT.NumberPoints2_1=numel(rho_t);

function tcol=init_mesh(t)
global OCMATCONT

arcposition=getarcposition(t);
stagenumber=uint8(length(arcposition));
arcnumber=uint8(zeros(1,stagenumber));
for ii=1:stagenumber
    arcnumber(ii)=size(arcposition{ii},2);
end
OCMATCONT.stagenumber=stagenumber;
OCMATCONT.arcnumber=arcnumber;

N=zeros(1,sum(arcnumber));
ctr=0;
for ii=1:stagenumber
    N(ctr+(1:arcnumber(ii)))=arcposition{ii}(2,:)-arcposition{ii}(1,:)+1;
    ctr=ctr+arcnumber(ii);
end
OCMATCONT.meshNumber=N;

tcol=makecollocationmesh(t);

function coeff=init_coefficients(sol)
global OCMATCONT
% to determine the Runge-Kutta coefficients the solution path is
% interpolated at the mesh with m+1 (number of collocation points)
% uniformly distributed grid points at each interval [t_{i-1} t_i] for
% first order and m uniformly distributed grid points at each
% interval [t_{i-1} t_i] for zero order components.
% Therefore the matrix A (given by the psi and psi0 values is evaluated at
% the unit interval with m+1 and m uniformly distributed grid points.
%
% The coeff

% componentnumber ... total number of first and zeror order components

coeff=zeros(sum(OCMATCONT.meshNumber-1)*((OCMATCONT.CollocationNumber+1)*sum(OCMATCONT.firstordercomponentnumber)+OCMATCONT.CollocationNumber*sum(OCMATCONT.zeroordercomponentnumber))+OCMATCONT.freeparameternum,1);
 %tfine=zeros(1,sum(OCMATCONT.meshNumber-1)*(OCMATCONT.CollocationNumber+1)+OCMATCONT.arcnumber);
ctr=0;
ctr1=0;
ctr2=0;
ctr3=0;
addidx=0;
tst=0;
firstorderidx=zeros(OCMATCONT.CollocationNumber+1,OCMATCONT.firstordercomponentnumber*OCMATCONT.diffarcposition(2,OCMATCONT.arcnumber));
zeroorderidx=zeros(OCMATCONT.CollocationNumber,OCMATCONT.zeroordercomponentnumber*OCMATCONT.diffarcposition(2,OCMATCONT.arcnumber));
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

fn=fieldnames(OCMATCONT);
l1=length(fn);
bcnum=OCMATCONT.arcnumber*OCMATCONT.firstordercomponentnumber+OCMATCONT.freeparameternum-1; % assumption that the number of components remain equal for the different arcs

addFpos=bcnum;
OCMATCONT.firstorderFinterioridx=zeros(OCMATCONT.firstordercomponentnumber*OCMATCONT.CollocationNumber*sum(OCMATCONT.meshNumber-1),1);
OCMATCONT.zeroorderFinterioridx=zeros(OCMATCONT.zeroordercomponentnumber*OCMATCONT.CollocationNumber*sum(OCMATCONT.meshNumber-1),1);
OCMATCONT.transitionFidx=zeros(OCMATCONT.firstordercomponentnumber*sum(OCMATCONT.meshNumber-2),1);

arcpart4firstorderFinterioridx=OCMATCONT.firstordercomponentnumber*OCMATCONT.CollocationNumber*cumsum(OCMATCONT.meshNumber-1);
OCMATCONT.arcpart4firstorderFinterioridx=[1 arcpart4firstorderFinterioridx(1:end-1)+1;arcpart4firstorderFinterioridx];

arcpart4zeroorderFinterioridx=OCMATCONT.zeroordercomponentnumber*OCMATCONT.CollocationNumber*cumsum(OCMATCONT.meshNumber-1);
OCMATCONT.arcpart4zeroorderFinterioridx=[1 arcpart4zeroorderFinterioridx(1:end-1)+1;arcpart4zeroorderFinterioridx];

arcpart4transitionFidx=OCMATCONT.firstordercomponentnumber*cumsum(OCMATCONT.meshNumber-2);
OCMATCONT.arcpart4transitionFidx=[1 arcpart4transitionFidx(1:end-1)+1;arcpart4transitionFidx];

for arc=1:OCMATCONT.arcnumber
    coeffarc=zeros(OCMATCONT.componentnumber,(OCMATCONT.CollocationNumber+1)*(OCMATCONT.meshNumber(arc)-1));
    coeffidx=zeros(OCMATCONT.componentnumber,(OCMATCONT.CollocationNumber+1)*(OCMATCONT.meshNumber(arc)-1));
    
    % equation index, total number of equations for first order ODEs and
    % AEs for each arc
    eqnum=OCMATCONT.componentnumber*OCMATCONT.CollocationNumber*(OCMATCONT.meshNumber(arc)-1)+OCMATCONT.firstordercomponentnumber*(OCMATCONT.meshNumber(arc)-2);
    Fidx=(1:eqnum);
    
    % equations at the mesh points
    transitionFidx=repmat([zeros(OCMATCONT.componentnumber*OCMATCONT.CollocationNumber,1);ones(OCMATCONT.firstordercomponentnumber,1)],OCMATCONT.meshNumber(arc)-1,1);
    transitionFidx=transitionFidx(1:eqnum);
    transitionFidx=find(transitionFidx==1);
    Fidx(transitionFidx)=[];
    Fidx=reshape(addFpos+Fidx,OCMATCONT.componentnumber,[]);
    OCMATCONT.transitionFidx(ctrtF+(1:OCMATCONT.firstordercomponentnumber*(OCMATCONT.meshNumber(arc)-2)))=addFpos+transitionFidx;
    OCMATCONT.firstorderFinterioridx(ctrfoFi+(1:OCMATCONT.firstordercomponentnumber*OCMATCONT.CollocationNumber*(OCMATCONT.meshNumber(arc)-1)))=Fidx(OCMATCONT.firstordercoordinate,:);
    OCMATCONT.zeroorderFinterioridx(ctrzoFi+(1:OCMATCONT.zeroordercomponentnumber*OCMATCONT.CollocationNumber*(OCMATCONT.meshNumber(arc)-1)))=Fidx(OCMATCONT.zeroordercoordinate,:);
    addFpos=addFpos+eqnum;
    ctrfoFi=ctrfoFi+OCMATCONT.firstordercomponentnumber*OCMATCONT.CollocationNumber*(OCMATCONT.meshNumber(arc)-1);
    ctrzoFi=ctrzoFi+OCMATCONT.zeroordercomponentnumber*OCMATCONT.CollocationNumber*(OCMATCONT.meshNumber(arc)-1);
    ctrtF=ctrtF+OCMATCONT.firstordercomponentnumber*(OCMATCONT.meshNumber(arc)-2);
    
    tst=tst+(OCMATCONT.CollocationNumber+1)*(OCMATCONT.meshNumber(arc)-1);
    t=sol.x(OCMATCONT.arcpositiondata(1,arc):OCMATCONT.arcpositiondata(2,arc));
    X=sol.y(:,OCMATCONT.arcpositiondata(1,arc):OCMATCONT.arcpositiondata(2,arc));
    tnew=sol.xnew(OCMATCONT.arcposition(1,arc):OCMATCONT.arcposition(2,arc));
    h=tnew(2:end)-tnew(1:end-1);

    psival=zeros(OCMATCONT.CollocationNumber,OCMATCONT.CollocationNumber+1);
    interp_t0=0:1/OCMATCONT.CollocationNumber:1;
    % evaluating first order collocation polynomials at interpolation times
    for ii=1:OCMATCONT.CollocationNumber
        psival(ii,1:OCMATCONT.CollocationNumber+1)=polyval(OCMATCONT.FirstOrderCollocationPolynomial(ii,:),interp_t0);
    end

    X1=X(OCMATCONT.firstordercoordinate,:);
    X0=X(OCMATCONT.zeroordercoordinate,:);

    interp_t1=[tnew(1:end-1);repmat(h/OCMATCONT.CollocationNumber,OCMATCONT.CollocationNumber-1,1)];
    interp_t1=cumsum(interp_t1);
    interp_t1=interp_t1(:).';
    interp_t1((OCMATCONT.meshNumber(arc)-1)*OCMATCONT.CollocationNumber+1)=tnew(end);
    
    if OCMATCONT.CollocationNumber>2
        interp_t0=[tnew(1:end-1);repmat(h/(OCMATCONT.CollocationNumber-1),OCMATCONT.CollocationNumber-2,1)];
        interp_t0=cumsum(interp_t0);
        interp_t0=interp_t0(:).';
        interp_t0((OCMATCONT.meshNumber(arc)-1)*(OCMATCONT.CollocationNumber-1)+1)=tnew(end);
    else
        % the interpolation is done at the mesh points
        interp_t0=tnew;
    end
    interp_X1=interp1(t,X1.',interp_t1,OCMATCONT.InterpolationMethod).';

    if size(X0,1)==1
        interp_X0=interp1(t,X0,interp_t0,OCMATCONT.InterpolationMethod);
    else
        interp_X0=interp1(t,X0.',interp_t0,OCMATCONT.InterpolationMethod).';
    end

    A=ones(OCMATCONT.CollocationNumber);
    for ii=1:OCMATCONT.CollocationNumber
        A(:,ii)=polyval(OCMATCONT.ZeroOrderCollocationPolynomial(ii,:),0:1/(OCMATCONT.CollocationNumber-1):1).';
    end
    locctr=0;
    locctr1=1;
    coeff0=NaN(OCMATCONT.zeroordercomponentnumber,(OCMATCONT.CollocationNumber+1)*(OCMATCONT.meshNumber(arc)-1));
    for ii=1:OCMATCONT.meshNumber(arc)-1
        coeff0(:,locctr1+(1:(OCMATCONT.CollocationNumber)))=(A\interp_X0(:,locctr+(1:OCMATCONT.CollocationNumber)).').';
        locctr1=locctr1+OCMATCONT.CollocationNumber+1;
        locctr=locctr+OCMATCONT.CollocationNumber-1;
    end

    Atmp=ones(OCMATCONT.CollocationNumber+1);
    locctr=0;
    locctr1=0;
    coeff1=zeros(OCMATCONT.firstordercomponentnumber,(OCMATCONT.CollocationNumber+1)*(OCMATCONT.meshNumber(arc)-1));
    for ii=1:OCMATCONT.meshNumber(arc)-1
        A=Atmp;
        A(:,2:OCMATCONT.CollocationNumber+1)=h(ii)*psival.';
        coeff1(:,locctr1+(1:(OCMATCONT.CollocationNumber+1)))=(A\interp_X1(:,locctr+(1:OCMATCONT.CollocationNumber+1)).').';
        locctr=locctr+OCMATCONT.CollocationNumber;
        locctr1=locctr1+OCMATCONT.CollocationNumber+1;
    end
    coeffarc(OCMATCONT.firstordercoordinate,:)=coeff1;
    coeffarc(OCMATCONT.zeroordercoordinate,:)=coeff0;
    coeffarc=coeffarc(:);
    idx=isnan(coeffarc);
    coeffarc(idx)=[];

    coeff(ctr+(1:(OCMATCONT.meshNumber(arc)-1)*((OCMATCONT.CollocationNumber+1)*OCMATCONT.firstordercomponentnumber+OCMATCONT.CollocationNumber*OCMATCONT.zeroordercomponentnumber)))=coeffarc;
    ctr=ctr+(OCMATCONT.meshNumber(arc)-1)*((OCMATCONT.CollocationNumber+1)*OCMATCONT.firstordercomponentnumber+OCMATCONT.CollocationNumber*OCMATCONT.zeroordercomponentnumber);

    tmptfine=[tnew(1:OCMATCONT.meshNumber(arc)-1);repmat(tnew(1:OCMATCONT.meshNumber(arc)-1),OCMATCONT.CollocationNumber,1)+repmat(h,OCMATCONT.CollocationNumber,1).*repmat(OCMATCONT.CollocationPoint(:),1,OCMATCONT.meshNumber(arc)-1)];
    tfine(ctr1+(1:(OCMATCONT.meshNumber(arc)-1)*(OCMATCONT.CollocationNumber+1)+1))=[tmptfine(:).' tnew(OCMATCONT.meshNumber(arc))];
    ctr1=ctr1+(OCMATCONT.meshNumber(arc)-1)*(OCMATCONT.CollocationNumber+1)+1;

    coeff0=ones(OCMATCONT.zeroordercomponentnumber,(OCMATCONT.CollocationNumber+1)*(OCMATCONT.meshNumber(arc)-1));
    coeff0(:,1:OCMATCONT.CollocationNumber+1:end)=0;
    coeff1=ones(OCMATCONT.firstordercomponentnumber,(OCMATCONT.CollocationNumber+1)*(OCMATCONT.meshNumber(arc)-1));

    coeffidx(OCMATCONT.firstordercoordinate,:)=coeff1;
    coeffidx(OCMATCONT.zeroordercoordinate,:)=coeff0;
    coeffidx=cumsum(coeffidx(:));
    coeffidx=reshape(coeffidx,OCMATCONT.componentnumber,[]);
    coeffidx(idx)=NaN;
    coeffidx=coeffidx+addidx;
    % index to retrieve first and zero order components from the coefficients
    firstorderidx(:,ctr2+(1:OCMATCONT.firstordercomponentnumber*(OCMATCONT.meshNumber(arc)-1)))=reshape(permute(reshape(coeffidx(OCMATCONT.firstordercoordinate,:),OCMATCONT.firstordercomponentnumber,OCMATCONT.CollocationNumber+1,OCMATCONT.meshNumber(arc)-1),[2 1 3]),OCMATCONT.CollocationNumber+1,[]);
    tmpzeroorderidx=reshape(permute(reshape(coeffidx(OCMATCONT.zeroordercoordinate,:),OCMATCONT.zeroordercomponentnumber,OCMATCONT.CollocationNumber+1,OCMATCONT.meshNumber(arc)-1),[2 1 3]),OCMATCONT.CollocationNumber+1,[]);
    tmpzeroorderidx(isnan(tmpzeroorderidx))=[];
    tmpzeroorderidx=reshape(tmpzeroorderidx,OCMATCONT.CollocationNumber,[]);
    zeroorderidx(:,ctr3+(1:OCMATCONT.zeroordercomponentnumber*(OCMATCONT.meshNumber(arc)-1)))=tmpzeroorderidx;
    addidx=addidx+max([firstorderidx(:);zeroorderidx(:)]);
    ctr2=ctr2+OCMATCONT.firstordercomponentnumber*(OCMATCONT.meshNumber(arc)-1);
    ctr3=ctr3+OCMATCONT.zeroordercomponentnumber*(OCMATCONT.meshNumber(arc)-1);
end
% 'firstorderidx' and 'zeororderidx' retrieve the coefficients of the order
% one and zero from the total coeffcient vector in the form (n ... number
% of first order components, m ... number of zero order components, N mesh
% size:
% firstorderidx: (tau_i T_{i,1} ... T_{i,m}) x n*(N-1) 
% zeroorderidx: (T_{i,1} ... T_{i,m}) x n*(N-1) 
OCMATCONT.firstorderidx=firstorderidx;
OCMATCONT.zeroorderidx=zeroorderidx;

% the arc partition for the 'firstorderidx' and 'zeororderidx'
arcpart4firstorderidx=cumsum(OCMATCONT.firstordercomponentnumber*(OCMATCONT.meshNumber-1));
OCMATCONT.arcpart4firstorderidx=[1 arcpart4firstorderidx(1:end-1)+1;arcpart4firstorderidx];
arcpart4zeroorderidx=cumsum(OCMATCONT.zeroordercomponentnumber*(OCMATCONT.meshNumber-1));
OCMATCONT.arcpart4zeroorderidx=[1 arcpart4zeroorderidx(1:end-1)+1;arcpart4zeroorderidx];

% the y and first order z coefficients are directly used for the system of
% equations. The according coefficients allow to retrieve these
% coefficients immediately from the total coefficient vector in the matrix
% form (k ... number collocation points): 
% ycoefficientidx : n x N-1
% firstorderzcoefficientidx : n x k*(N-1)
OCMATCONT.ycoefficientidx=reshape(OCMATCONT.firstorderidx(1,:),OCMATCONT.firstordercomponentnumber,[]); % index to retrieve the y coefficients
OCMATCONT.firstorderzcoefficientidx=reshape(permute(reshape(OCMATCONT.firstorderidx(2:OCMATCONT.CollocationNumber+1,:),OCMATCONT.CollocationNumber,OCMATCONT.firstordercomponentnumber,sum(OCMATCONT.meshNumber-1)),[2 1 3]),OCMATCONT.firstordercomponentnumber,[]); % index to retrieve the z coefficients for first order components

% the arc partition for the 'ycoefficientidx' and 'firstorderzcoefficientidx'
arcpart4ycoefficientidx=cumsum(OCMATCONT.meshNumber(arc)-1);
OCMATCONT.arcpart4ycoefficientidx=[1 arcpart4ycoefficientidx(1:end-1)+1;arcpart4ycoefficientidx];

arcpart4firstorderzcoefficientidx=cumsum(OCMATCONT.CollocationNumber*(OCMATCONT.meshNumber-1));
OCMATCONT.arcpart4firstorderzcoefficientidx=[1 arcpart4firstorderzcoefficientidx(1:end-1)+1;arcpart4firstorderzcoefficientidx];

coeff(ctr+(1:OCMATCONT.freeparameternum))=sol.parameters;

OCMATCONT.freeparameterindex=ctr+(1:OCMATCONT.freeparameternum); % position of the free parameters in the coefficient vector
OCMATCONT.continuationindex=ctr+OCMATCONT.freeparameternum;

fn=fieldnames(OCMATCONT);
l2=length(fn);

for ii=(l1+1):l2
    OCMATCONT.TMPDATA.(fn{ii})=OCMATCONT.(fn{ii});
end
OCMATCONT.bcidx=1:bcnum;
OCMATCONT.bcnum=bcnum;

function init_jacobian(sol)
global OCMATCONT

h=diff(sol.xnew);
% JBlock ... block of psi entries for first order components derived from
% the derivative with respect to c_ikj
% BasicJBlock ... ones derived from the derivative with respect to y_ik,
% JBlockidx ... is position index of JBlock in BasicJBlock, i.e.
% BasicJBlock(JBlockidx)=JBlock 
% 
% BasicJBlocknonzeroidx ... indices of nonzero elements of BasicJBlock
% BasicJBlockAddOneidx ... indices where a one is added to the derivative
% of F_r(Til) with respect to c_iks, since r=k and l=s.
OCMATCONT.JBlock=zeros(OCMATCONT.componentnumber*OCMATCONT.CollocationNumber,OCMATCONT.firstordercomponentnumber*OCMATCONT.CollocationNumber);
OCMATCONT.JBlockidx=OCMATCONT.firstordercomponentnumber+OCMATCONT.firstordercoordinate(:,ones(1,OCMATCONT.CollocationNumber))+OCMATCONT.componentnumber*repmat(0:OCMATCONT.CollocationNumber-1,OCMATCONT.firstordercomponentnumber,1);
OCMATCONT.JBlockidx=OCMATCONT.JBlockidx(:).';
OCMATCONT.BasicJBlock=zeros(OCMATCONT.componentnumber*OCMATCONT.CollocationNumber,OCMATCONT.CollocationNumber*OCMATCONT.componentnumber+OCMATCONT.firstordercomponentnumber);
OCMATCONT.BasicJBlock(:,1:OCMATCONT.firstordercomponentnumber)=1;

OCMATCONT.BasicJBlocknonzeroidx=zeros(1,OCMATCONT.zeroordercomponentnumber*OCMATCONT.componentnumber*(OCMATCONT.CollocationNumber-1));
OCMATCONT.BasicJBlockAddOneidx=zeros(1,OCMATCONT.firstordercomponentnumber*OCMATCONT.CollocationNumber);
OCMATCONT.BasicJBlockSize=[OCMATCONT.CollocationNumber*OCMATCONT.componentnumber,OCMATCONT.CollocationNumber*OCMATCONT.componentnumber+OCMATCONT.firstordercomponentnumber];
OCMATCONT.BasicJBlocknonzeroidx=[];
OCMATCONT.CBlock=zeros(OCMATCONT.firstordercomponentnumber,OCMATCONT.firstordercomponentnumber*OCMATCONT.CollocationNumber);
OCMATCONT.BasicCBlock=zeros(OCMATCONT.firstordercomponentnumber,OCMATCONT.CollocationNumber*OCMATCONT.componentnumber+2*OCMATCONT.firstordercomponentnumber);

OCMATCONT.BasicCBlock(:,1:OCMATCONT.firstordercomponentnumber)=eye(OCMATCONT.firstordercomponentnumber);
OCMATCONT.BasicCBlock(:,OCMATCONT.CollocationNumber*OCMATCONT.componentnumber+OCMATCONT.firstordercomponentnumber+(1:OCMATCONT.firstordercomponentnumber))=-eye(OCMATCONT.firstordercomponentnumber);
OCMATCONT.BasicCBlockSize=[OCMATCONT.firstordercomponentnumber,OCMATCONT.CollocationNumber*OCMATCONT.componentnumber+2*OCMATCONT.firstordercomponentnumber];
OCMATCONT.BasicCBlocknonzeroidx=[];

fn=fieldnames(OCMATCONT);
l1=length(fn);


OCMATCONT.BCbBlock=ones(OCMATCONT.firstordercomponentnumber+OCMATCONT.freeparameternum-1,OCMATCONT.firstordercomponentnumber+OCMATCONT.componentnumber*OCMATCONT.CollocationNumber);
OCMATCONT.BCbBlockidx=zeros(1,OCMATCONT.firstordercomponentnumber+OCMATCONT.componentnumber*OCMATCONT.CollocationNumber);

%tic
psival=OCMATCONT.FirstOrderCollocationValues(2:OCMATCONT.CollocationNumber+1,2:OCMATCONT.CollocationNumber+1);
%psival=psival.';
psival=psival(:).';
psival=psival(ones(OCMATCONT.componentnumber*OCMATCONT.firstordercomponentnumber,1),:);

ctr=0;
rows=1:OCMATCONT.componentnumber;
cols=1:OCMATCONT.firstordercomponentnumber;
tmp=zeros(OCMATCONT.componentnumber,OCMATCONT.firstordercomponentnumber);
tmp1=ones(OCMATCONT.componentnumber,OCMATCONT.zeroordercomponentnumber);
tmp2=cumsum([0 (OCMATCONT.firstordercoordinate(2:OCMATCONT.firstordercomponentnumber)-OCMATCONT.firstordercoordinate(1:OCMATCONT.firstordercomponentnumber-1)).']*OCMATCONT.componentnumber*OCMATCONT.CollocationNumber);

tmpb1=zeros(OCMATCONT.firstordercomponentnumber);
tmpb0=zeros(OCMATCONT.firstordercomponentnumber,OCMATCONT.zeroordercomponentnumber);
for ii=0:OCMATCONT.CollocationNumber-1
    for jj=0:OCMATCONT.CollocationNumber-1
        ctr=ctr+1;
        tmp(:)=psival(:,ctr);
        OCMATCONT.JBlock(ii*OCMATCONT.componentnumber+rows,jj*OCMATCONT.firstordercomponentnumber+cols)=tmp;
        if ii==jj
            OCMATCONT.BasicJBlock(ii*OCMATCONT.componentnumber+rows,OCMATCONT.firstordercomponentnumber+jj*OCMATCONT.componentnumber+OCMATCONT.zeroordercoordinate)=tmp1;
        end
    end
    OCMATCONT.BasicJBlockAddOneidx(ii*OCMATCONT.firstordercomponentnumber+1:(ii+1)*OCMATCONT.firstordercomponentnumber)=(OCMATCONT.firstordercomponentnumber+ii*OCMATCONT.componentnumber)*OCMATCONT.componentnumber*OCMATCONT.CollocationNumber+ii*OCMATCONT.componentnumber+OCMATCONT.firstordercoordinate(:).'+tmp2;
    
    OCMATCONT.CBlock(:,ii*OCMATCONT.firstordercomponentnumber+cols)=eye(OCMATCONT.firstordercomponentnumber)*OCMATCONT.FirstOrderCollocationValues(ii+2,OCMATCONT.CollocationNumber+2);
    
    tmpb1(:)=h(end-1)*OCMATCONT.FirstOrderCollocationValues(ii+2,OCMATCONT.CollocationNumber+2);
    tmpb0(:)=OCMATCONT.ZeroOrderCollocationValues(ii+1,OCMATCONT.CollocationNumber+2);
    OCMATCONT.BCbBlock(:,OCMATCONT.firstordercomponentnumber+ii*OCMATCONT.componentnumber+OCMATCONT.firstordercoordinate)=tmpb1;
    OCMATCONT.BCbBlock(:,OCMATCONT.firstordercomponentnumber+ii*OCMATCONT.componentnumber+OCMATCONT.zeroordercoordinate)=tmpb0;
end

BasicJBlock=OCMATCONT.BasicJBlock;
BasicCBlock=OCMATCONT.BasicCBlock;
BasicJBlock(:,OCMATCONT.JBlockidx)=OCMATCONT.JBlock;
BasicJBlock(OCMATCONT.BasicJBlockAddOneidx)=1;
BasicCBlock(:,OCMATCONT.JBlockidx)=OCMATCONT.CBlock;
OCMATCONT.BasicJBlocknonzeroidx=find(BasicJBlock).';
%OCMATCONT.BasicJBlockzeroidx=find(BasicJBlock==0).';
OCMATCONT.BasicCBlocknonzeroidx=find(BasicCBlock).';
%OCMATCONT.BasicCBlockzeroidx=find(BasicCBlock==0).';

OCMATCONT.BCbBlockidx=(OCMATCONT.meshNumber-2)*(OCMATCONT.componentnumber*OCMATCONT.CollocationNumber+OCMATCONT.firstordercomponentnumber)+ ...
    (1:OCMATCONT.firstordercomponentnumber+OCMATCONT.componentnumber*OCMATCONT.CollocationNumber);

OCMATCONT.JacobianNonzeroEntriesNumber=(OCMATCONT.meshNumber-1)*OCMATCONT.componentnumber*OCMATCONT.CollocationNumber*(OCMATCONT.firstordercomponentnumber*OCMATCONT.CollocationNumber+1)+ ...
    (OCMATCONT.meshNumber-2)*OCMATCONT.firstordercomponentnumber*(OCMATCONT.CollocationNumber+2)+ ...
    (OCMATCONT.firstordercomponentnumber+OCMATCONT.freeparameternum-1)*(OCMATCONT.componentnumber*OCMATCONT.CollocationNumber+ ...
    2*OCMATCONT.firstordercomponentnumber)+((OCMATCONT.meshNumber-1)*OCMATCONT.componentnumber*OCMATCONT.CollocationNumber+OCMATCONT.firstordercomponentnumber+ ...
    OCMATCONT.freeparameternum-1)*OCMATCONT.freeparameternum;

fn=fieldnames(OCMATCONT);
l2=length(fn);

for ii=(l1+1):l2
    OCMATCONT.TMPDATA.(fn{ii})=OCMATCONT.(fn{ii});
end



