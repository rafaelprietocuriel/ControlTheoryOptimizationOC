function res=calcres_test(tmesh,coeff,dae,bc)
global OCMATCONT OCMATFTE OCBVP
tic
[tmesh,x,y,z1,freepar,modelpar]=rearr(tmesh,coeff);
arcnumber=OCBVP.arcnumber;
arcpositioncollocationmesh=OCBVP.arcpositioncollocationmesh;

res=zeros(OCMATCONT.continuationindex,1);
res(OCBVP.bcidx) = bc(x(:,OCBVP.arcpositioncollocationmesh(1,:)),x(:,OCBVP.arcpositioncollocationmesh(2,:)),freepar,modelpar);

for arc = 1:arcnumber
    % retrieves indices of the collocation points (idxmesh) and at the
    % meshpoints (idxbd)
    idxmesh=arcpositioncollocationmesh(1,arc):arcpositioncollocationmesh(2,arc)-1;
    idxbd=idxmesh(1:OCBVP.CollocationNumber+1:end);
    idxmesh(idxbd)=[];
    
    
    % z_ijk=f() at the collocation points the first order components
    % satisfy the ODE
    fz=dae(tmesh(idxmesh),x(:,idxmesh),freepar,modelpar,arc);
    res(OCBVP.firstorderFinterioridx(OCBVP.arcpart4firstorderFinterioridx(1,arc):OCBVP.arcpart4firstorderFinterioridx(2,arc)))= ...
        fz(OCBVP.firstordercoordinate,:)-z1(:,OCBVP.arcpart4firstorderzcoefficientidx(1,arc):OCBVP.arcpart4firstorderzcoefficientidx(2,arc));
    
    %y=P(0) at the interior mesh points the first order components are
    %continuous
    res(OCBVP.transitionFidx(OCBVP.arcpart4transitionFidx(1,arc):OCBVP.arcpart4transitionFidx(2,arc)))= ...
        y(:,(OCBVP.arcpart4ycoefficientidx(1,arc)+1):OCBVP.arcpart4ycoefficientidx(2,arc))-x(OCBVP.firstordercoordinate,idxbd(2:end));
    
    % f()=0 at the collocation points the zero order equations are
    % satisfied
    res(OCBVP.zeroorderFinterioridx(OCBVP.arcpart4zeroorderFinterioridx(1,arc):OCBVP.arcpart4zeroorderFinterioridx(2,arc)))= ...
        fz(OCBVP.zeroordercoordinate,:);
end
tst=toc;
disp(tst)


order=zeros(OCBVP.componentnumber,1);
order(OCBVP.firstordercoordinate)=1;
psi1=zeros(OCBVP.CollocationNumber,max(order)+OCBVP.CollocationNumber,max(order));
for ord=1:max(order)
    for i=1:OCBVP.CollocationNumber
        psi1(i,1+max(order)-ord:OCBVP.CollocationNumber+max(order),ord)=Psi(i,OCBVP.CollocationPoint,ord);
    end
end
psival1=zeros(max(order),OCBVP.CollocationNumber,OCBVP.CollocationNumber+2);
for ord=1:max(order)
    for i=1:OCBVP.CollocationNumber
        %evaluation of psi at the collocation points,0 and 1
        psival1(ord,i,1:OCBVP.CollocationNumber)=polyval(psi1(i,:,ord),OCBVP.CollocationPoint(1:OCBVP.CollocationNumber));
        psival1(ord,i,OCBVP.CollocationNumber+1)=polyval(psi1(i,:,ord),1);
        psival1(ord,i,OCBVP.CollocationNumber+2)=polyval(psi1(i,:,ord),0);
    end
end
tic
fval=functionFDF('F','bvps2_nonlincapacc_probdef',coeff,OCBVP.t,psival1,psi1,OCBVP.CollocationPoint,[]);
tst=toc;
disp(tst)
max(abs(fval))


% ------------------------------------------------------
function [tmesh,x,y,z1,freepar,modelpar]=rearr(tmesh,coeff)
global OCMATCONT OCMATFTE OCBVP

modelpar=OCMATCONT.modelparameter;
freepar=coeff(OCMATCONT.continuationindex); % parameter values
if ~isempty(OCMATCONT.freeparameter)
    modelpar(OCMATCONT.freeparameter)=freepar(OCMATCONT.freeparametercoordinate);
end
y=coeff(OCBVP.ycoefficientidx); % y coefficients
z1=coeff(OCBVP.firstorderzcoefficientidx); % z coefficients for the first order components
x=coeff2points(coeff,'collocationgrid'); % returns the values of the states approximated by the collocation polynomials at the mesh with collocation points


function prod=Psi(n,rho,nr)
% calculates the Lagrangepolynomials and its integrals
global OCBVP

prod=1;
for ii=1:OCBVP.CollocationNumber
    if (ii~=n)
        prod=conv(prod,[1 -rho(ii)])/(rho(n)-rho(ii));
    end
end
if nr==1
    prod=polyint(prod);
end
