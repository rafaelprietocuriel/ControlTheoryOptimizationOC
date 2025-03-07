function sol=adapt4daesolver(sol,order,collMethod,collocnum,interpMethod)
% 
% for the computation the solution components are reordered such that the
% first order components are at the beginnig and thereafter the zero order
% components. The sol.y components are in the initial ordering.

global OCMATCONT

t=sol.x;
X=sol.y;
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

if isfield(sol,'arcposition')
    arcposition=sol.arcposition;
else
    arcposition=[1;length(t)];
end
arcnumber=size(arcposition,2);
OCMATCONT.arcposition=arcposition;
OCMATCONT.arcnumber=arcnumber;
% in the actual implementation it is assumed that each arc consists of the
% same components, therefore variables like 'order',
% 'firstordercoordinate', etc. are arc-independent

if isempty(t) || isempty(X)
    return
end

if length(t)~=size(X,2)
    return
end

if isfield(sol,'data') && isfield(sol.data,'coeff') &&  isfield(sol.data,'tfine')
    coeff0=sol.data.coeff;
    tfine0=sol.data.tfine;
else
    coeff0=[];
end

if nargin<5
    opt=defaultocoptions;
    interpMethod=getocoptions(opt,'SBVPOC','InterpolationMethod');
end

if nargin<4
    collocnum=getocoptions(opt,'SBVPOC','CollocationNumber');
end

if nargin<3
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

% [ordernew,orderidx]=sort(order,'descend'); % order(orderindex)=ordernew
% tmp=1:length(order);
% [dum,reorderidx]=sort(tmp(orderidx)); % ordernew(reorderindex)=order
OCMATCONT.freeparameternum=length(freepar);

OCMATCONT.firstordercoordinate=find(order==1);
OCMATCONT.zeroordercoordinate=find(order==0);
OCMATCONT.componentnumber=length(order);
OCMATCONT.firstordercomponentnumber=length(OCMATCONT.firstordercoordinate);
OCMATCONT.zeroordercomponentnumber=length(OCMATCONT.zeroordercoordinate);
%OCMATCONT.reorderidx=reorderidx;
%OCMATCONT.orderidx=orderidx;

OCMATCONT.CollocationNumber=collocnum;
OCMATCONT.CollocationMethod=collMethod;
OCMATCONT.InterpolationMethod=interpMethod;

N=arcposition(2,:)-arcposition(1,:)+1;
OCMATCONT.meshNumber=N;
idx=cumsum((N-1)*(collocnum+1)+1);
OCMATCONT.arcpositioncollocationmesh=[1 idx(1:arcnumber-1)+1;idx];

h=t(2:end)-t(1:end-1);
cumsumNm1=cumsum(N-1);
h(cumsumNm1(1:end-1)+(1:arcnumber-1))=[]; % remove time difference between two adjacent arcs 

tmp=arcposition(2,:)-(1:arcnumber);
OCMATCONT.diffarcposition=[1 tmp(1:end-1)+1;tmp];

OCMATCONT.t=t;
OCMATCONT.h=h;

H=repmat(OCMATCONT.h,OCMATCONT.firstordercomponentnumber,1);
H=repmat(H(:).',OCMATCONT.CollocationNumber,1);
OCMATCONT.H=H;

rho=getcollocationpoints(collMethod,collocnum);
psi=zeros(collocnum,1+collocnum);
for ii=1:collocnum
    psi(ii,1:collocnum+1)=Psi(ii,rho,1);
end

psival=zeros(collocnum+1,collocnum+2);
psival(1,:)=1; % first row corresponds to values of the polynomial derived by the integration of psi.   
for ii=1:collocnum
    %evaluation of psi at the collocation points,0 and 1 for order 1
    %components
    %psival(ii+1,1)=polyval(psi(ii,:),0);
    psival(ii+1,1+(1:collocnum))=polyval(psi(ii,:),rho(1:collocnum));
    psival(ii+1,collocnum+2)=polyval(psi(ii,:),1);
end

%evaluation of psi at the collocation points,0 and 1 for order 0
%components
psi0=zeros(collocnum,collocnum);
for ii=1:collocnum
    psi0(ii,:)=Psi(ii,rho,0);
end
psi0val=zeros(collocnum,collocnum+2);
for ii=1:collocnum
    %evaluation of psi at the collocation points, 0 and 1 for order 0
    %components
    psi0val(ii,1)=polyval(psi0(ii,:),0);
    psi0val(ii,1+(1:collocnum))=polyval(psi0(ii,:),rho(1:collocnum));
    psi0val(ii,collocnum+2)=polyval(psi0(ii,:),1);
end

OCMATCONT.CollocationPoint=rho;
OCMATCONT.FirstOrderCollocationPolynomial=psi;
OCMATCONT.ZeroOrderCollocationPolynomial=psi0;
OCMATCONT.FirstOrderCollocationValues=psival;
OCMATCONT.ZeroOrderCollocationValues=psi0val;

% calculate Runge-Kutta coefficients from solution data
%tic
% coeff=zeros(sum(N-1)*((collocnum+1)*OCMATCONT.firstordercomponentnumber+collocnum*OCMATCONT.zeroordercomponentnumber)+OCMATCONT.freeparameternum,1);
% tfine=zeros(1,sum(N-1)*(collocnum+1)+arcnumber);
% ctr=0;
% ctr1=0;
% for ii=1:arcnumber
%     coeff(ctr+(1:(N(ii)-1)*((collocnum+1)*OCMATCONT.firstordercomponentnumber+collocnum*OCMATCONT.zeroordercomponentnumber)))=init_coefficients(t(arcposition(1,ii):arcposition(2,ii)),X(:,arcposition(1,ii):arcposition(2,ii)),ii);
%     ctr=ctr+(N(ii)-1)*((collocnum+1)*OCMATCONT.firstordercomponentnumber+collocnum*OCMATCONT.zeroordercomponentnumber);
%     tmptfine=[t(arcposition(1,ii):arcposition(2,ii)-1);repmat(t(arcposition(1,ii):arcposition(2,ii)-1),collocnum,1)+repmat(h(OCMATCONT.diffarcposition(1,ii):OCMATCONT.diffarcposition(2,ii)),collocnum,1).*repmat(rho(:),1,N(ii)-1)];
%     tfine(ctr1+(1:(N(ii)-1)*(collocnum+1)+1))=[tmptfine(:).' t(arcposition(2,ii))];
%     ctr1=ctr1+(N(ii)-1)*(collocnum+1)+1;
% end
[coeff,tfine]=init_coefficients(sol);
%tst=toc;
%disp(tst)
% if all(size(coeff)==size(coeff0)) && length(tfine0)==length(tfine) && max(abs(tfine-tfine0))<1e-12
%     % if the collocation mesh of the solution sol coincides with the
%     % computed mesh and the size of the Runge-Kutta coefficients is equal,
%     % the original data is used.
%     coeff=coeff0;
%     tfine=tfine0;
% end
sol.data.coeff=coeff;
sol.data.tfine=tfine;

OCMATCONT.tfine=sol.data.tfine;
% tic
% initProfile.initialMesh=t;
% initProfile.initialValues=X;
% initProfile=initial_coefficients(['bvps2_nonlincapacc_probdef'],t,initProfile,rho,1);
% YY=coeffToValues(initProfile.initialCoeff,initProfile.initialMesh,order,rho);
% sol.data.initt=t;
% sol.data.YY=YY;

% tst=toc;
% disp(tst)
%

%Y=coeff2points(coeff,'grid');
Yf=coeff2points(coeff,'collocationgrid');
sol.data.Yf=Yf;

Y=coeff2points(coeff,'grid');
sol.data.Y=Y;

% number of boundary conditions
OCMATCONT.nBCs=OCMATCONT.firstordercomponentnumber+length(freepar);

init_jacobian();

function prod=Psi(n,rho,nr)
% calculates the Lagrangepolynomials and its integrals
global OCMATCONT

prod=1;
for ii=1:OCMATCONT.CollocationNumber
    if (ii~=n)
        prod=conv(prod,[1 -rho(ii)])/(rho(n)-rho(ii));
    end
end
if nr==1
    prod=polyint(prod);
end


function rho=getcollocationpoints(type,k)
if (strcmp(type,'gauss'))
    %psi(:,:,i) for n-th order psi(:,:,n-i) equals i-th derivative
    %for example for 3-rd order:
    %psi(:,:,1) ... second derivative
    %psi(:,:,2) ... first derivative
    %psi(:,:,3) ... no derivative

    switch k
        case 1
            rho=[0.5];
            psi(1,:,1)=[0 1 0];
            psi(1,:,2)=[0.5 0 0];
        case 2
            %C: rho: ci in Gaussian collocation, zeros of Legendre Polynomials shiftet to [0,1]
            %C: transformation from [-1,1] to [0,1]: y= 1/2x+1/2
            rho=[0.211324865405187117745425609748  0.788675134594812882254574390252];
            psi(1,:,1)=[0 -0.866025403784438646763723170755  1.36602540378443864676372317076 0];
            psi(2,:,1)=[0  0.866025403784438646763723170755   -0.366025403784438646763723170754 0];
            %0.th derivative
            %firs coll.point
            psi(1,:,2)=[-0.288675134594812882254574390252  0.683012701892219323381861585378 0 0];
            %second coll.point
            psi(2,:,2)=[0.288675134594812882254574390252  -0.183012701892219323381861585377 0 0];
        case 3
            rho=[1/2-1/10*15^(1/2) 1/2 1/2+1/10*15^(1/2)];
            psi(1,:,1)=[0 10/9 -5/3-1/6*15^(1/2) 5/6+1/6*15^(1/2) 0];
            psi(2,:,1)=[0 -20/9 10/3 -2/3 0];
            psi(3,:,1)=[0 10/9 -5/3+1/6*15^(1/2) 5/6-1/6*15^(1/2) 0];
            psi(1,:,2)=[5/18 -5/9-1/18*15^(1/2) 5/12+1/12*15^(1/2) 0 0];
            psi(2,:,2)=[-5/9 10/9 -1/3 0 0];
            psi(3,:,2)=[5/18 -5/9+1/18*15^(1/2) 5/12-1/12*15^(1/2) 0 0];
        case 4

            rho=[1/2-1/70*(525+70*30^(1/2))^(1/2) 1/2-1/70*(525-70*30^(1/2))^(1/2) 1/2+1/70*(525-70*30^(1/2))^(1/2) 1/2+1/70*(525+70*30^(1/2))^(1/2)];
            psi(1,:,1)=[0 -245/24/(525+70*30^(1/2))^(1/2)*30^(1/2) 7/36*(105+(525+70*30^(1/2))^(1/2))/(525+70*30^(1/2))^(1/2)*30^(1/2) ...
                -7/24*(45+(525+70*30^(1/2))^(1/2)+30^(1/2))/(525+70*30^(1/2))^(1/2)*30^(1/2) ...
                1/120*(35+(525+70*30^(1/2))^(1/2))*(10+30^(1/2))*30^(1/2)/(525+70*30^(1/2))^(1/2) 0];
            psi(2,:,1)=[0 245/24/(525-70*30^(1/2))^(1/2)*30^(1/2) -7/36*(105+(525-70*30^(1/2))^(1/2))/(525-70*30^(1/2))^(1/2)*30^(1/2) ...
                -7/24*(-45+30^(1/2)-(525-70*30^(1/2))^(1/2))/(525-70*30^(1/2))^(1/2)*30^(1/2) ...
                1/120*(35+(525-70*30^(1/2))^(1/2))*(-10+30^(1/2))*30^(1/2)/(525-70*30^(1/2))^(1/2) 0];
            psi(3,:,1)=[0 -245/24/(525-70*30^(1/2))^(1/2)*30^(1/2) -7/36*(-105+(525-70*30^(1/2))^(1/2))/(525-70*30^(1/2))^(1/2)*30^(1/2) ...
                7/24*(-45+30^(1/2)+(525-70*30^(1/2))^(1/2))/(525-70*30^(1/2))^(1/2)*30^(1/2) ...
                1/120*(-35+(525-70*30^(1/2))^(1/2))*(-10+30^(1/2))*30^(1/2)/(525-70*30^(1/2))^(1/2) 0];
            psi(4,:,1)=[0 245/24/(525+70*30^(1/2))^(1/2)*30^(1/2) 7/36*(-105+(525+70*30^(1/2))^(1/2))/(525+70*30^(1/2))^(1/2)*30^(1/2) ...
                7/24*(45-(525+70*30^(1/2))^(1/2)+30^(1/2))/(525+70*30^(1/2))^(1/2)*30^(1/2) ...
                1/120*(-35+(525+70*30^(1/2))^(1/2))*(10+30^(1/2))*30^(1/2)/(525+70*30^(1/2))^(1/2) 0];
            psi(1,:,2)=[-49/24/(525+70*30^(1/2))^(1/2)*30^(1/2) 7/144*(105+(525+70*30^(1/2))^(1/2))/(525+70*30^(1/2))^(1/2)*30^(1/2) ...
                -7/72*(45+(525+70*30^(1/2))^(1/2)+30^(1/2))/(525+70*30^(1/2))^(1/2)*30^(1/2) ...
                1/240*(35+(525+70*30^(1/2))^(1/2))*(10+30^(1/2))*30^(1/2)/(525+70*30^(1/2))^(1/2) 0 0];
            psi(2,:,2)=[49/24/(525-70*30^(1/2))^(1/2)*30^(1/2) -7/144*(105+(525-70*30^(1/2))^(1/2))/(525-70*30^(1/2))^(1/2)*30^(1/2) ...
                7/72*(45-30^(1/2)+(525-70*30^(1/2))^(1/2))*30^(1/2)/(525-70*30^(1/2))^(1/2) ...
                1/240*(35+(525-70*30^(1/2))^(1/2))*(-10+30^(1/2))*30^(1/2)/(525-70*30^(1/2))^(1/2) 0 0];
            psi(3,:,2)=[-49/24/(525-70*30^(1/2))^(1/2)*30^(1/2) -7/144*(-105+(525-70*30^(1/2))^(1/2))/(525-70*30^(1/2))^(1/2)*30^(1/2) ...
                7/72*(-45+30^(1/2)+(525-70*30^(1/2))^(1/2))/(525-70*30^(1/2))^(1/2)*30^(1/2) ...
                1/240*(-35+(525-70*30^(1/2))^(1/2))*(-10+30^(1/2))*30^(1/2)/(525-70*30^(1/2))^(1/2) 0 0];
            psi(4,:,2)=[49/24/(525+70*30^(1/2))^(1/2)*30^(1/2) 7/144*(-105+(525+70*30^(1/2))^(1/2))/(525+70*30^(1/2))^(1/2)*30^(1/2) ...
                7/72*(45-(525+70*30^(1/2))^(1/2)+30^(1/2))/(525+70*30^(1/2))^(1/2)*30^(1/2) ...
                1/240*(-35+(525+70*30^(1/2))^(1/2))*(10+30^(1/2))*30^(1/2)/(525+70*30^(1/2))^(1/2) 0 0];

        case 5

            rho=[1/2-1/42*(245+14*70^(1/2))^(1/2) 1/2-1/42*(245-14*70^(1/2))^(1/2) 1/2 1/2+1/42*(245-14*70^(1/2))^(1/2) 1/2+1/42*(245+14*70^(1/2))^(1/2)];
            psi(1,:,1)=[0 567/25/(35+2*70^(1/2))*70^(1/2) -27/40*(84+(245+14*70^(1/2))^(1/2))/(35+2*70^(1/2))*70^(1/2) ...
                3/20*(343+9*(245+14*70^(1/2))^(1/2)+2*70^(1/2))/(35+2*70^(1/2))*70^(1/2) ...
                -3/280*(1911+77*(245+14*70^(1/2))^(1/2)+42*70^(1/2)+(245+14*70^(1/2))^(1/2)*70^(1/2))/(35+2*70^(1/2))*70^(1/2) ...
                3/280*(21+(245+14*70^(1/2))^(1/2))*(14+70^(1/2))*70^(1/2)/(35+2*70^(1/2)) 0];
            psi(2,:,1)=[0 567/25/(-35+2*70^(1/2))*70^(1/2) -27/40*(84+(245-14*70^(1/2))^(1/2))/(-35+2*70^(1/2))*70^(1/2) ...
                -3/20*(-343+2*70^(1/2)-9*(245-14*70^(1/2))^(1/2))/(-35+2*70^(1/2))*70^(1/2) ...
                3/280*(-1911+42*70^(1/2)-77*(245-14*70^(1/2))^(1/2)+(245-14*70^(1/2))^(1/2)*70^(1/2))/(-35+2*70^(1/2))*70^(1/2) ...
                -3/280*(21+(245-14*70^(1/2))^(1/2))*(-14+70^(1/2))*70^(1/2)/(-35+2*70^(1/2)) 0];
            psi(3,:,1)=[0 336/25 -168/5 1232/45 -112/15 8/15 0];
            psi(4,:,1)=[0 567/25/(-35+2*70^(1/2))*70^(1/2) 27/40*(-84+(245-14*70^(1/2))^(1/2))/(-35+2*70^(1/2))*70^(1/2) ...
                -3/20*(-343+2*70^(1/2)+9*(245-14*70^(1/2))^(1/2))/(-35+2*70^(1/2))*70^(1/2) ...
                -3/280*(1911-42*70^(1/2)-77*(245-14*70^(1/2))^(1/2)+(245-14*70^(1/2))^(1/2)*70^(1/2))/(-35+2*70^(1/2))*70^(1/2) ...
                3/280*(-21+(245-14*70^(1/2))^(1/2))*(-14+70^(1/2))*70^(1/2)/(-35+2*70^(1/2)) 0];
            psi(5,:,1)=[0 567/25/(35+2*70^(1/2))*70^(1/2) 27/40*(-84+(245+14*70^(1/2))^(1/2))/(35+2*70^(1/2))*70^(1/2) ...
                -3/20*(-343+9*(245+14*70^(1/2))^(1/2)-2*70^(1/2))/(35+2*70^(1/2))*70^(1/2) ...
                3/280*(-1911+77*(245+14*70^(1/2))^(1/2)-42*70^(1/2)+(245+14*70^(1/2))^(1/2)*70^(1/2))/(35+2*70^(1/2))*70^(1/2) ...
                -3/280*(-21+(245+14*70^(1/2))^(1/2))*(14+70^(1/2))*70^(1/2)/(35+2*70^(1/2)) 0];
            psi(1,:,2)=[189/50/(35+2*70^(1/2))*70^(1/2) -27/200*(84+(245+14*70^(1/2))^(1/2))/(35+2*70^(1/2))*70^(1/2) ...
                3/80*(343+9*(245+14*70^(1/2))^(1/2)+2*70^(1/2))/(35+2*70^(1/2))*70^(1/2) ...
                -1/280*(1911+77*(245+14*70^(1/2))^(1/2)+42*70^(1/2)+(245+14*70^(1/2))^(1/2)*70^(1/2))/(35+2*70^(1/2))*70^(1/2) ...
                3/560*(21+(245+14*70^(1/2))^(1/2))*(14+70^(1/2))*70^(1/2)/(35+2*70^(1/2)) 0 0];
            psi(2,:,2)=[189/50/(-35+2*70^(1/2))*70^(1/2) -27/200*(84+(245-14*70^(1/2))^(1/2))/(-35+2*70^(1/2))*70^(1/2) ...
                3/80*(343-2*70^(1/2)+9*(245-14*70^(1/2))^(1/2))/(-35+2*70^(1/2))*70^(1/2) ...
                1/280*(-1911+42*70^(1/2)-77*(245-14*70^(1/2))^(1/2)+(245-14*70^(1/2))^(1/2)*70^(1/2))/(-35+2*70^(1/2))*70^(1/2) ...
                -3/560*(21+(245-14*70^(1/2))^(1/2))*(-14+70^(1/2))*70^(1/2)/(-35+2*70^(1/2)) 0 0];
            psi(3,:,2)=[56/25 -168/25 308/45 -112/45 4/15 0 0];
            psi(4,:,2)=[189/50/(-35+2*70^(1/2))*70^(1/2) 27/200*(-84+(245-14*70^(1/2))^(1/2))/(-35+2*70^(1/2))*70^(1/2) -3/80*(-343+2*70^(1/2)+9*(245-14*70^(1/2))^(1/2))/(-35+2*70^(1/2))*70^(1/2) ...
                -1/280*(1911-42*70^(1/2)-77*(245-14*70^(1/2))^(1/2)+(245-14*70^(1/2))^(1/2)*70^(1/2))/(-35+2*70^(1/2))*70^(1/2) ...
                3/560*(-21+(245-14*70^(1/2))^(1/2))*(-14+70^(1/2))*70^(1/2)/(-35+2*70^(1/2)) 0 0];
            psi(5,:,2)=[189/50/(35+2*70^(1/2))*70^(1/2) 27/200*(-84+(245+14*70^(1/2))^(1/2))/(35+2*70^(1/2))*70^(1/2) ...
                3/80*(343-9*(245+14*70^(1/2))^(1/2)+2*70^(1/2))*70^(1/2)/(35+2*70^(1/2)) ...
                1/280*(-1911+77*(245+14*70^(1/2))^(1/2)-42*70^(1/2)+(245+14*70^(1/2))^(1/2)*70^(1/2))/(35+2*70^(1/2))*70^(1/2) ...
                -3/560*(-21+(245+14*70^(1/2))^(1/2))*(14+70^(1/2))*70^(1/2)/(35+2*70^(1/2)) 0 0];
        case 6

            rho=[0.33765242898423986094e-1 0.16939530676686774317 0.38069040695840154568 0.61930959304159845432 0.83060469323313225683 0.96623475710157601391];
            psi(1,:,1)=[0 -8.1412617290086756177 28.978672208695686824 -40.408362140839295291 27.785390555068868256 -9.6944498478780709320 1.5656732001510719330 0];
            psi(2,:,1)=[0 24.533738740695531559 -83.334379226361945114 107.81105264377978693 -64.863346973441684340 16.973778445028729194 -.94046284317634892902 0];
            psi(3,:,1)=[0 -36.168340495472103373 113.68329746902200059 -130.85406519770904404 65.101388388983683720 -12.145253252968680079 .61693005543048870860 0];
            psi(4,:,1)=[0 36.168340495472103373 -103.32674550381061965 104.96268528468059170 -44.848707621074554040 7.6576120141334378913 -.37922770211461375460 0];
            psi(5,:,1)=[0 -24.533738740695531559 63.868053217811244240 -59.145237622403034740 23.711846151968643422 -3.9123422341959200109 .19180001403866795482 0];
            psi(6,:,1)=[0 8.1412617290086756177 -19.868898165356366883 17.633927032490995438 -6.8865705015049570236 1.1206548758805039360 -.54712724329265912892e-1 0];
            psi(1,:,2)=[-1.1630373898583822311 4.8297787014492811373 -8.0816724281678590582 6.9463476387672170640 -3.2314832826260236440 .78283660007553596650 0 0];
            psi(2,:,2)=[3.5048198200993616513 -13.889063204393657519 21.562210528755957386 -16.215836743360421085 5.6579261483429097313 -.47023142158817446451 0 0];
            psi(3,:,2)=[-5.1669057850674433390 18.947216244837000098 -26.170813039541808808 16.275347097245920930 -4.0484177509895600263 .30846502771524435430 0 0];
            psi(4,:,2)=[5.1669057850674433390 -17.221124250635103275 20.992537056936118340 -11.212176905268638510 2.5525373380444792971 -.18961385105730687730 0 0];
            psi(5,:,2)=[-3.5048198200993616513 10.644675536301874040 -11.829047524480606948 5.9279615379921608555 -1.3041140780653066703 .95900007019333977410e-1 0 0];
            psi(6,:,2)=[1.1630373898583822311 -3.3114830275593944805 3.5267854064981990876 -1.7216426253762392559 .37355162529350131200 -.27356362164632956446e-1 0 0];
        case 7
            rho=[.254460438286207377369051579761e-1, .129234407200302780068067613360, .297077424311301416546696793962, .500000000000000000000000000000, .702922575688698583453303206038, .870765592799697219931932386640, .974553956171379262263094842024];
        case 8
            rho=[.198550717512318841582195657153e-1, .101666761293186630204223031762, .237233795041835507091130475405, .408282678752175097530261928820, .591717321247824902469738071180, .762766204958164492908869524595, .898333238706813369795776968238, .980144928248768115841780434285];
        case 9
            rho=[.159198802461869550822118985482e-1, .819844463366821028502851059651e-1, .193314283649704801345648980329, .337873288298095535480730992678, .500000000000000000000000000000, .662126711701904464519269007322, .806685716350295198654351019671, .918015553663317897149714894035, .984080119753813044917788101452];
        case 10
            rho=[.130467357414141399610179939578e-1, .674683166555077446339516557883e-1, .160295215850487796882836317443, .283302302935376404600367028417, .425562830509184394557586999435, .574437169490815605442413000565, .716697697064623595399632971583, .839704784149512203117163682557, .932531683344492255366048344212, .986953264258585860038982006042];
        case 11
            rho=[.108856709269715035980309994386e-1, .564687001159523504624211153480e-1, .134923997212975337953291873984, .240451935396594092037137165271, .365228422023827513834234007300, .500000000000000000000000000000, .634771577976172486165765992700, .759548064603405907962862834729, .865076002787024662046708126016, .943531299884047649537578884652, .989114329073028496401969000561];
        case 12
            rho=[.921968287664037465472545492536e-2, .479413718147625716607670669405e-1, .115048662902847656481553083394, .206341022856691276351648790530, .316084250500909903123654231678, .437383295744265542263779315268, .562616704255734457736220684732, .683915749499090096876345768322, .793658977143308723648351209470, .884951337097152343518446916606, .952058628185237428339232933060, .990780317123359625345274545075];
        case 13
            rho=[.790847264070592526358527559645e-2, .412008003885110173967260817496e-1, .992109546333450436028967552086e-1, .178825330279829889678007696502, .275753624481776573561043573936, .384770842022432602967235939451, .500000000000000000000000000000, .615229157977567397032764060549, .724246375518223426438956426064, .821174669720170110321992303498, .900789045366654956397103244791, .958799199611488982603273918250, .992091527359294074736414724404];
        case 14
            rho=[.685809565159383057920136664797e-2, .357825581682132413318044303111e-1, .863993424651175034051026286748e-1, .156353547594157264925990098490, .242375681820922954017354640724, .340443815536055119782164087916, .445972525646328168966877674890, .554027474353671831033122325110, .659556184463944880217835912084, .757624318179077045982645359276, .843646452405842735074009901510, .913600657534882496594897371325, .964217441831786758668195569689, .993141904348406169420798633352];
        case 15
            rho=[.600374098975728575521714070669e-2, .313633037996470478461205261449e-1, .758967082947863918996758396129e-1, .137791134319914976291906972693, .214513913695730576231386631373, .302924326461218315051396314509, .399402953001282738849685848303, .500000000000000000000000000000, .600597046998717261150314151697, .697075673538781684948603685491, .785486086304269423768613368627, .862208865680085023708093027307, .924103291705213608100324160387, .968636696200352952153879473855, .993996259010242714244782859293];
    end
end
if(strcmp(type,'lobatto'))
    if k<2 || k>15
        err('equations_err4');
    end
    switch k
        case 2
            rho=[0,1];
        case 3
            rho=[0., .500000000000000000000000000000, 1.];
        case 4
            rho=[0., .276393202250021030359082633127, .723606797749978969640917366873, 1.];
        case 5
            rho=[0., .172673164646011428100853771876, .500000000000000000000000000000, .827326835353988571899146228124, 1.];
        case 6
            rho=[0., .117472338035267653574498513019, .357384241759677451842924502980, .642615758240322548157075497020, .882527661964732346425501486981, 1.];
        case 7
            rho=[0., .84888051860716535063983893017e-1, .265575603264642893098114059046, .500000000000000000000000000000, .734424396735357106901885940954, .915111948139283464936016106983, 1.];
        case 8
            rho=[0., .641299257451966923312771193897e-1, .204149909283428848927744634301, .395350391048760565615671369827, .604649608951239434384328630173, .795850090716571151072255365699, .935870074254803307668722880610, 1.00000000000000000000000000000];
        case 9
            rho=[0., .501210022942699213438273777908e-1, .161406860244631123277057286454, .318441268086910920644623965646, .500000000000000000000000000000, .681558731913089079355376034354, .838593139755368876722942713546, .949878997705730078656172622209, 1.00000000000000000000000000000];
        case 10
            rho=[0., .402330459167705930855336695888e-1, .130613067447247462498446912570, .261037525094777752169412453634, .417360521166806487686890117021, .582639478833193512313109882979, .738962474905222247830587546366, .869386932552752537501553087430, .959766954083229406914466330411, 1.00000000000000000000000000000];
        case 11
            rho=[0., .329992847959704328338629319503e-1, .107758263168427790688791091946, .217382336501897496764518015261, .352120932206530304284044242220, .500000000000000000000000000000, .647879067793469695715955757780, .782617663498102503235481984739, .892241736831572209311208908054, .967000715204029567166137068050, 1.00000000000000000000000000000];
        case 12
            rho=[0., .275503638885588882962099308484e-1, .903603391779966608256792091415e-1, .183561923484069661168797572778, .300234529517325533867825104217, .431723533572536222567969072130, .568276466427463777432030927870, .699765470482674466132174895783, .816438076515930338831202427222, .909639660822003339174320790858, .972449636111441111703790069152, 1.00000000000000000000000000000];
        case 13
            rho=[0., .233450766789180440515472676223e-1, .768262176740638415670371964506e-1, .156905765459121286963620480217, .258545089454331899126531383182, .375356534946880003715663149813, .500000000000000000000000000000, .624643465053119996284336850187, .741454910545668100873468616818, .843094234540878713036379519783, .923173782325936158432962803549, .976654923321081955948452732378, 1.00000000000000000000000000000];
        case 14
            rho=[0., .200324773663695493224499189923e-1, .660994730848263744998898985459e-1, .135565700454336929707663799740, .224680298535676472341688647070, .328637993328643577478048298179, .441834065558148066170611645132, .558165934441851933829388354868, .671362006671356422521951701821, .775319701464323527658311352930, .864434299545663070292336200260, .933900526915173625500110101454, .979967522633630450677550081008, 1.00000000000000000000000000000];
        case 15
            rho=[0., .173770367480807136020743039652e-1, .574589778885118505872991842589e-1, .118240155024092399647940762012, .196873397265077144438235030682, .289680972643163759539051530631, .392323022318102880887160276864, .500000000000000000000000000000, .607676977681897119112839723136, .710319027356836240460948469369, .803126602734922855561764969318, .881759844975907600352059237988, .942541022111488149412700815741, .982622963251919286397925696035, 1.00000000000000000000000000000];
    end
end
if(strcmp(type,'uniform'))
    rho = 1/(k+1):1/(k+1):k/(k+1);
end

function [coeff,tfine]=init_coefficients(sol)
global OCMATCONT
% to determine the Runge-Kutta coefficients the solution path is
% interpolated at the mesh with m+1 (number of collocation points)
% uniformly distributed grid points at each interval [t_{i-1} t_i] for
% first order  and m uniformly distributed grid points at each
% interval [t_{i-1} t_i] for zero order components.
% Therefore the matrix A (given by the psi and psi0 values is evaluated at
% the unit interval with m+1 and m uniformly distributed grid points.
%
% The coeff

% componentnumber ... total number of first and zeror order components

coeff=zeros(sum(OCMATCONT.meshNumber-1)*((OCMATCONT.CollocationNumber+1)*OCMATCONT.firstordercomponentnumber+OCMATCONT.CollocationNumber*OCMATCONT.zeroordercomponentnumber)+OCMATCONT.freeparameternum,1);
tfine=zeros(1,sum(OCMATCONT.meshNumber-1)*(OCMATCONT.CollocationNumber+1)+OCMATCONT.arcnumber);
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

OCMATCONT.bcidx=1:bcnum;
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
    t=sol.x(OCMATCONT.arcposition(1,arc):OCMATCONT.arcposition(2,arc));
    X=sol.y(:,OCMATCONT.arcposition(1,arc):OCMATCONT.arcposition(2,arc));
    h=t(2:end)-t(1:end-1);

    psival=zeros(OCMATCONT.CollocationNumber,OCMATCONT.CollocationNumber+1);
    interp_t0=0:1/OCMATCONT.CollocationNumber:1;
    % evaluating first order collocation polynomials at interpolation times
    for ii=1:OCMATCONT.CollocationNumber
        psival(ii,1:OCMATCONT.CollocationNumber+1)=polyval(OCMATCONT.FirstOrderCollocationPolynomial(ii,:),interp_t0);
    end

    X1=X(OCMATCONT.firstordercoordinate,:);
    X0=X(OCMATCONT.zeroordercoordinate,:);

    interp_t1=[t(1:end-1);repmat(h/OCMATCONT.CollocationNumber,OCMATCONT.CollocationNumber-1,1)];
    interp_t1=cumsum(interp_t1);
    interp_t1=interp_t1(:).';
    interp_t1((OCMATCONT.meshNumber(arc)-1)*OCMATCONT.CollocationNumber+1)=t(end);
    
    if OCMATCONT.CollocationNumber>2
        interp_t0=[t(1:end-1);repmat(h/(OCMATCONT.CollocationNumber-1),OCMATCONT.CollocationNumber-2,1)];
        interp_t0=cumsum(interp_t0);
        interp_t0=interp_t0(:).';
        interp_t0((OCMATCONT.meshNumber(arc)-1)*(OCMATCONT.CollocationNumber-1)+1)=t(end);
    else
        % the interpolation is done at the mesh points
        interp_t0=t;
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

    tmptfine=[t(1:OCMATCONT.meshNumber(arc)-1);repmat(t(1:OCMATCONT.meshNumber(arc)-1),OCMATCONT.CollocationNumber,1)+repmat(h,OCMATCONT.CollocationNumber,1).*repmat(OCMATCONT.CollocationPoint(:),1,OCMATCONT.meshNumber(arc)-1)];
    tfine(ctr1+(1:(OCMATCONT.meshNumber(arc)-1)*(OCMATCONT.CollocationNumber+1)+1))=[tmptfine(:).' t(OCMATCONT.meshNumber(arc))];
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
arcpart4firstorderidx=cumsum(OCMATCONT.firstordercomponentnumber*(OCMATCONT.meshNumber(arc)-1));
OCMATCONT.arcpart4firstorderidx=[1 arcpart4firstorderidx(1:end-1)+1;arcpart4firstorderidx];
arcpart4zeroorderidx=cumsum(OCMATCONT.zeroordercomponentnumber*(OCMATCONT.meshNumber(arc)-1));
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

function init_jacobian()
global OCMATCONT

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

OCMATCONT.CBlock=zeros(OCMATCONT.firstordercomponentnumber,OCMATCONT.firstordercomponentnumber*OCMATCONT.CollocationNumber);
OCMATCONT.BasicCBlock=zeros(OCMATCONT.firstordercomponentnumber,OCMATCONT.CollocationNumber*OCMATCONT.componentnumber+2*OCMATCONT.firstordercomponentnumber);

OCMATCONT.BasicCBlock(:,1:OCMATCONT.firstordercomponentnumber)=eye(OCMATCONT.firstordercomponentnumber);
OCMATCONT.BasicCBlock(:,OCMATCONT.CollocationNumber*OCMATCONT.componentnumber+OCMATCONT.firstordercomponentnumber+(1:OCMATCONT.firstordercomponentnumber))=-eye(OCMATCONT.firstordercomponentnumber);
OCMATCONT.BasicCBlockSize=[OCMATCONT.firstordercomponentnumber,OCMATCONT.CollocationNumber*OCMATCONT.componentnumber+2*OCMATCONT.firstordercomponentnumber];

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
    
    tmpb1(:)=OCMATCONT.h(OCMATCONT.meshNumber-1)*OCMATCONT.FirstOrderCollocationValues(ii+2,OCMATCONT.CollocationNumber+2);
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

