function out=evaluteDF_LC(coeff,tangent)
limitcycleData=limitcycleargument();
ocMFVar=occontargument();

c=feval(ocMFVar.bvpfile,'c');%boundary values inside the interval
[y,z,p,mode]=rearrcoefflc(coeff);
[yref,zref]=rearrcoefflc(limitcycleData.coeffref);
switch mode
    case 1
        gridsize=limitcycleData.gridsize;
        tau=limitcycleData.tau;
        N=limitcycleData.N;
        ROW=limitcycleData.ROW;
        COL=limitcycleData.COL;
        dimsparse=limitcycleData.dimsparse;
        nonzeroent=limitcycleData.nonzeroent;

    case 2
        gridsize=limitcycleData.gridsize2;
        tau=limitcycleData.tau2;
        N=limitcycleData.N2;
        ROW=limitcycleData.ROW2;
        COL=limitcycleData.COL2;
        dimsparse=limitcycleData.dimsparse2;
        nonzeroent=limitcycleData.nonzeroent2;
end
%JACOBIAN

% Dimensioning of non-zero entries
% Boundary conditions + Matrix
% gen_c ... General count for sparse matrix efficiency

% OUTPUT=zeros(nonzeroent,1);
% unitmatrix=eye(ocMFVar.n);
% fz=zeros(ocMFVar.n,2,ocMFVar.ncol);
% zerostemplate=zeros(1,ocMFVar.n);
% onestemplate=ones(1,ocMFVar.n);
% helptemplate=ones(1,ocMFVar.ncol);
% unitmatrix2=zeros(ocMFVar.n);
% unitmatrix2(:,:,2)=unitmatrix;
OUTPUT=zeros(nonzeroent,1);
unitmatrix=limitcycleData.template.eyen;
fz=limitcycleData.template.zerosn2ncol;
zerostemplate=limitcycleData.template.zerosn;
onestemplate=limitcycleData.template.onesn;
helptemplate=limitcycleData.template.onesncol;
unitmatrix2=limitcycleData.template.zerosnn;
unitmatrix2(:,:,2)=unitmatrix;

DCtot_=[ones(N-1,1) kron(gridsize(1:N-1).',limitcycleData.psival(1,:,ocMFVar.ncolp1)) -ones(N-1,1) ];

if (N==1)
    error('Number of time intervals less than two.')
end
if ~isempty(c)
    error('Only TPBVP are considered at this place.')
end
%----Start Boundary conditions
%R_0

fzb=Poptimiert0(N,gridsize,y,z,ocMFVar,limitcycleData);
DRb_=feval(ocMFVar.bvpfile,'DR','b1',ocMFVar.nppcoord,[],[],y(:,1,1),fzb,[],[],p).';
DI=evaluteDIC(y,z,zref,limitcycleData,ocMFVar,N,gridsize);
DR=[feval(ocMFVar.bvpfile,'DR','a1',ocMFVar.nppcoord,[],[],y(:,1,1),fzb,[],[],p).'; ...
    DRb_; ...
    kron(gridsize(N)*limitcycleData.psival(:,:,ocMFVar.ncolp1),DRb_(:,:).').'; ...
    feval(ocMFVar.bvpfile,'DpR',[],[],[],[],y(:,1,1),fzb,[],[],p).'];
gen_c=(ocMFVar.npp-1)*(ocMFVar.n*ocMFVar.ncolp2+ocMFVar.numparameter);
OUTPUT(1:gen_c)=DR(:);
gen_c_start=gen_c+1;
gen_c=gen_c+N*(ocMFVar.n*ocMFVar.ncol+ocMFVar.n);
OUTPUT(gen_c_start:gen_c)=DI(1:end-1).';
%----------End Boundary Conditions

%----------Start J_i with C_i (notation according to manual)

for i=1:N
    fz(:,1,:)=Poptimiert2(i,gridsize,y,z,ocMFVar,limitcycleData);
    fz(:,2,:)=z(:,:,i);
    fztest=fz(:,1,:);
    fztest=fztest(:,:);
    Dgvect_=feval(ocMFVar.bvpfile,'Dgvect',1,[],0,fztest(:,ocMFVar.ncolcoord),[],[],[],[],p);
    for k=ocMFVar.ncolcoord
        template=helptemplate;
        template(k)=2;
%         Dg_=feval(ocMFVar.bvpfile,'Dg',1,[],tau(i,k),fz(:,:,k),[],[],[],[],p);
        %Dg_=Dgtest_(:,limitcycleData.jacpatterng(k,:));
        DG=Dgvect_(:,limitcycleData.jacpatterng(k,:),ones(1,ocMFVar.ncol));%Dg_(:,:,ones(1,ocMFVar.ncol));
        %DG=Dg_(:,:,ones(1,ocMFVar.ncol));
        DG=gridsize(i)*DG.*limitcycleData.psival2(k+zerostemplate,onestemplate,:)+unitmatrix2(:,:,template);
        Dpg=feval(ocMFVar.bvpfile,'Dpg',[],[],tau(i,k),fz(:,:,k),[],[],[],[],p);
        DG=[Dgvect_(:,limitcycleData.jacpatterng(k,:)) DG(:,:) Dpg].';
        %DG=[Dg_ DG(:,:) Dpg].';
        gen_c_start=gen_c+1;
        gen_c=gen_c+ocMFVar.n*(ocMFVar.n*ocMFVar.ncolp1+ocMFVar.numparameter);
        OUTPUT(gen_c_start:gen_c)=DG(:);
    end
    if i<N
        % derivative of continuity conditions
        gen_c_start=gen_c+1;
        gen_c=gen_c+ocMFVar.n*ocMFVar.ncolp2;
        DC_=DCtot_(i,:).';
        DC_=DC_(:,onestemplate);
        OUTPUT(gen_c_start:gen_c)=DC_(:);
    end
end

out=sparse(ROW,COL,OUTPUT,dimsparse,dimsparse);

%--------End J_i with C_i

%Change of the derivation of the ocMFVar.n+1 th row when using
%pathfollowing

if ~isempty(tangent)
    if size(tangent,1)==2
        % replace row ocMFVar.n+1 by the previously detected
        % tangent vector
        out(ocMFVar.n+1,:)=tangent(2,:);
    end
    if size(tangent,1)==3
        out(ocMFVar.n+1,N*(ocMFVar.n*ocMFVar.ncol+ocMFVar.n)+1)=1;
    end
end
%--------Ende Jacobian

function ret=Poptimiert0(i,h,y,z,ocMFVar,limitcycleData)
%Polynomial opitmized version

ret=y(:,1,i)+sum(h(i)*z(:,:,i).*limitcycleData.psivalm(:,:,ocMFVar.ncolp1),2);

function ret=Poptimiert2(i,h,y,z,ocMFVar,limitcycleData)
%Polynomial opitmized version

zidx=i+zeros(1,ocMFVar.ncol);
ret=y(:,1,zidx)+sum(h((i))*z(:,:,zidx).*limitcycleData.psivalm(:,:,ocMFVar.ncolcoord),2);
ret=ret(:,:);
