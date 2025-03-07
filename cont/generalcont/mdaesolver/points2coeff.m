function coeff=points2coeff(sol,tnew)
global OCMATCONT
% returns the coefficients of the solution sol, for the meshgrid tnew

idx=find(diff(sol.x)==0);
arcpositiondata=[1 idx+1;idx length(sol.x)];
idx=find(diff(tnew)==0);
arcposition=[1 idx+1;idx length(tnew)];
arcnumber=size(arcposition,2);

meshnumber=diff(arcposition)+1;
coeff=zeros(sum(meshnumber-1)*((OCMATCONT.CollocationNumber+1)*OCMATCONT.firstordercomponentnumber+OCMATCONT.CollocationNumber*OCMATCONT.zeroordercomponentnumber)+OCMATCONT.freeparameternum,1);

ctr=0;
for arc=1:arcnumber
    t=sol.x(arcpositiondata(1,arc):arcpositiondata(2,arc));
    X=sol.y(:,arcpositiondata(1,arc):arcpositiondata(2,arc));
    tnewarc=tnew(arcposition(1,arc):arcposition(2,arc));
    h=tnewarc(2:end)-tnewarc(1:end-1);

    psival=zeros(OCMATCONT.CollocationNumber,OCMATCONT.CollocationNumber+1);
    interp_t0=0:1/OCMATCONT.CollocationNumber:1;
    % evaluating first order collocation polynomials at interpolation times
    for ii=1:OCMATCONT.CollocationNumber
        psival(ii,1:OCMATCONT.CollocationNumber+1)=polyval(OCMATCONT.FirstOrderCollocationPolynomial(ii,:),interp_t0);
    end

    X1=X(OCMATCONT.firstordercoordinate,:);
    X0=X(OCMATCONT.zeroordercoordinate,:);

    interp_t1=[tnewarc(1:end-1);repmat(h/OCMATCONT.CollocationNumber,OCMATCONT.CollocationNumber-1,1)];
    interp_t1=cumsum(interp_t1);
    interp_t1=interp_t1(:).';
    interp_t1((meshnumber(arc)-1)*OCMATCONT.CollocationNumber+1)=tnewarc(end);
    
    if OCMATCONT.CollocationNumber>2
        interp_t0=[tnewarc(1:end-1);repmat(h/(OCMATCONT.CollocationNumber-1),OCMATCONT.CollocationNumber-2,1)];
        interp_t0=cumsum(interp_t0);
        interp_t0=interp_t0(:).';
        interp_t0((meshnumber(arc)-1)*(OCMATCONT.CollocationNumber-1)+1)=tnewarc(end);
    else
        % the interpolation is done at the mesh points
        interp_t0=tnewarc;
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
    coeff0=NaN(OCMATCONT.zeroordercomponentnumber,(OCMATCONT.CollocationNumber+1)*(meshnumber(arc)-1));
    for ii=1:meshnumber(arc)-1
        coeff0(:,locctr1+(1:(OCMATCONT.CollocationNumber)))=(A\interp_X0(:,locctr+(1:OCMATCONT.CollocationNumber)).').';
        %coeff0(:,locctr1+(1:(OCMATCONT.CollocationNumber)))=(A\interp_X0(:,locctr+(1:OCMATCONT.CollocationNumber)).');
        locctr1=locctr1+OCMATCONT.CollocationNumber+1;
        locctr=locctr+OCMATCONT.CollocationNumber-1;
    end

    Atmp=ones(OCMATCONT.CollocationNumber+1);
    locctr=0;
    locctr1=0;
    coeff1=zeros(OCMATCONT.firstordercomponentnumber,(OCMATCONT.CollocationNumber+1)*(meshnumber(arc)-1));
    for ii=1:meshnumber(arc)-1
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

    coeff(ctr+(1:(meshnumber(arc)-1)*((OCMATCONT.CollocationNumber+1)*OCMATCONT.firstordercomponentnumber+OCMATCONT.CollocationNumber*OCMATCONT.zeroordercomponentnumber)))=coeffarc;
    ctr=ctr+(meshnumber(arc)-1)*((OCMATCONT.CollocationNumber+1)*OCMATCONT.firstordercomponentnumber+OCMATCONT.CollocationNumber*OCMATCONT.zeroordercomponentnumber);
end