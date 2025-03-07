function output=evaluteF_LC(coeff,tangent)
limitcycleData=limitcycleargument();
ocMFVar=occontargument();

%----------------------------------------------
%Call data from bvpfile

[y,z,p,mode]=rearrcoefflc(coeff);
[yref,zref]=rearrcoefflc(limitcycleData.coeffref);
switch mode
    case 1
        gridsize=limitcycleData.gridsize;
        tau=limitcycleData.tau;
        N=limitcycleData.N;

    case 2
        gridsize=limitcycleData.gridsize2;
        tau=limitcycleData.tau2;
        N=limitcycleData.N2;
end
%in the rows of x1(1:limitcycleData.N).'*ones(1,ocMFVar.ncol) is xi
% in the columns are xi, ...xN

output=zeros(length(coeff),1);
jj=0;
Rx=feval(ocMFVar.bvpfile,'R',[],[],[],[],y(:,1,1),Poptimiert0(N,gridsize,y,z,ocMFVar,limitcycleData),[],[],p,[],[],[],[],[],ocMFVar);
if ocMFVar.numparameter>0
    output(1+jj:ocMFVar.npp-1+jj,1)=Rx;
    jj=jj+ocMFVar.npp;
end
output(jj)=evaluteIC(y,z,zref,ocMFVar,limitcycleData,N,gridsize);
if ~isempty(tangent)
    if size(tangent,1)==2
        %output(sum(ordnung)+1,1)=sum((a.'-tangent(1,:)).*tangent(2,:)); % maybe this can be set to zero, see MATCONT newtcorr
        output(ocMFVar.np1,1)=0; % set to zero, see MATCONT newtcorr
    end
    if size(tangent,1)==3
        output(ocMFVar.np1,1)=p(1)-tangent(3,1);
    end
end

%-----------------------------------------------
%Evaluation of the equations for the solver
fz=zeros(ocMFVar.n,2,ocMFVar.ncol);
for i=1:N
    fz(:,1,:)=Poptimiert2(i,gridsize,y,z,ocMFVar,limitcycleData); % limitcycleData.psi wird immer von der Ordnung der entsprechenden Gleichung gewählt
    fz(:,2,:)=z(:,:,i);
    for k=ocMFVar.ncolcoord
        % Dynamik muiss an den Kollokationspunkten erfüllt
        % sein
        g_=feval(ocMFVar.bvpfile,'g',[],[],tau(i,k),fz(:,:,k),[],[],[],[],p,[],[],[],[],[],ocMFVar);
        output(jj+ocMFVar.coord,1)=g_;
        jj=jj+ocMFVar.n;
    end
    if i<N
        %u=Poptimiert0(i,1:ocMFVar.n,x1((i+1)+1),ocMFVar.ncol,x1,gridsize,y,z,limitcycleData.psival,ocMFVar.ncol+1)-y(1:ocMFVar.n,1,i+2);
        output(jj+ocMFVar.coord)=Poptimiert0(i,gridsize,y,z,ocMFVar,limitcycleData)-y(:,1,i+1);
        jj=jj+ocMFVar.n;
    end
end

% function ret=Poptimiert(abl,i,komp,x,ocMFVar.ncol,x1,gridsize,y,z,psival,ord,k)
% %Polynomial opitmized version
% faktorielle=[1 1 2 6 24 120 720 5040 40320 362880 3628800 39916800 479001600 6227020800];
% 
% ret=sum(((x-x1((i)+1)).^([abl+1:ord]-abl-1))./(faktorielle([abl+1:ord]-abl)).*y(komp,abl+1:ord,(i)+1))+...
%     sum(gridsize((i)+1)^(ord-abl).*(z(komp,1:ocMFVar.ncol,(i)+1).*psival(ord-(abl+1)+1,1:ocMFVar.ncol,k)));

function ret=Poptimiert0(i,h,y,z,ocMFVar,limitcycleData)
%Polynomial opitmized version

ret=y(:,1,i)+sum(h(i)*z(:,:,i).*limitcycleData.psivalm(:,:,ocMFVar.ncolp1),2);

function ret=Poptimiert2(i,h,y,z,ocMFVar,limitcycleData)
%Polynomial opitmized version
%limitcycleData=manifoldargument;
%ocMFVar=occontargument();

zidx=i+zeros(1,ocMFVar.ncol);
ret=y(:,1,zidx)+sum(h((i))*z(:,:,zidx).*limitcycleData.psivalm(:,:,ocMFVar.ncolcoord),2);
ret=ret(:,:);
