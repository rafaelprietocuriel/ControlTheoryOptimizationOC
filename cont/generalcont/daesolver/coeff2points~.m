function X=coeff2points(coeff,flag)
global OCMATCONT

if nargin==1
    % grid with collocation points
    flag='collocationgrid';
end

switch flag
    case 'collocationgrid'
        % the components are computed als the polynomial values evaluated
        % at rho1 ... rho_m 1

        coeff1=coeff(OCMATCONT.firstorderidx);
        coeff1(2:OCMATCONT.CollocationNumber+1,:)=OCMATCONT.H.*coeff1(2:OCMATCONT.CollocationNumber+1,:);
        coeff0=coeff(OCMATCONT.zeroorderidx);
        
        X=zeros(OCMATCONT.componentnumber,sum(OCMATCONT.meshNumber-1)*(OCMATCONT.CollocationNumber+1)+OCMATCONT.arcnumber);
        ctr1=0;
        ctr0=0;
        ctr2=0;
        for ii=1:OCMATCONT.arcnumber
            arc_coeff1=coeff1(:,ctr1+(1:((OCMATCONT.meshNumber(ii)-1)*OCMATCONT.firstordercomponentnumber)));
            arc_coeff0=coeff0(:,ctr0+(1:((OCMATCONT.meshNumber(ii)-1)*OCMATCONT.zeroordercomponentnumber)));
            ctr1=ctr1+(OCMATCONT.meshNumber(ii)-1)*OCMATCONT.firstordercomponentnumber;
            ctr0=ctr0+(OCMATCONT.meshNumber(ii)-1)*OCMATCONT.zeroordercomponentnumber;
            X1=reshape(permute(reshape(OCMATCONT.FirstOrderCollocationValues(:,2:OCMATCONT.CollocationNumber+2).'*arc_coeff1,OCMATCONT.CollocationNumber+1,OCMATCONT.firstordercomponentnumber,[]),[2 1 3]),OCMATCONT.firstordercomponentnumber,[]);
            X0=reshape(permute(reshape(OCMATCONT.ZeroOrderCollocationValues(:,2:OCMATCONT.CollocationNumber+2).'*arc_coeff0,OCMATCONT.CollocationNumber+1,OCMATCONT.zeroordercomponentnumber,[]),[2 1 3]),OCMATCONT.zeroordercomponentnumber,[]);
            X10=OCMATCONT.FirstOrderCollocationValues(:,1).'*arc_coeff1(:,1:OCMATCONT.firstordercomponentnumber);
            X00=OCMATCONT.ZeroOrderCollocationValues(:,1).'*arc_coeff0(:,1:OCMATCONT.zeroordercomponentnumber);
            X(OCMATCONT.firstordercoordinate,ctr2+(1:((OCMATCONT.meshNumber(ii)-1)*(OCMATCONT.CollocationNumber+1)+1)))=[X10(:) X1];
            X(OCMATCONT.zeroordercoordinate,ctr2+(1:((OCMATCONT.meshNumber(ii)-1)*(OCMATCONT.CollocationNumber+1)+1)))=[X00(:) X0];
            ctr2=ctr2+((OCMATCONT.meshNumber(ii)-1)*(OCMATCONT.CollocationNumber+1)+1);
        end
    case 'collocationgrid_'
        % the components are computed als the polynomial values evaluated
        % at 0 rho1 ... rho_m
        X=zeros(OCMATCONT.componentnumber,sum(OCMATCONT.meshNumber-1)*(OCMATCONT.CollocationNumber+1)+OCMATCONT.arcnumber);
        ctr1=0;
        ctr0=0;
        ctr2=0;
        for ii=1:OCMATCONT.arcnumber
            arc_coeff1=coeff1(:,ctr1+(1:((OCMATCONT.meshNumber(ii)-1)*OCMATCONT.firstordercomponentnumber)));
            arc_coeff0=coeff0(:,ctr0+(1:((OCMATCONT.meshNumber(ii)-1)*OCMATCONT.zeroordercomponentnumber)));
            ctr1=ctr1+(OCMATCONT.meshNumber(ii)-1)*OCMATCONT.firstordercomponentnumber;
            ctr0=ctr0+(OCMATCONT.meshNumber(ii)-1)*OCMATCONT.zeroordercomponentnumber;
            X1=reshape(permute(reshape(OCMATCONT.FirstOrderCollocationValues(:,1:OCMATCONT.CollocationNumber+1).'*arc_coeff1,OCMATCONT.CollocationNumber+1,OCMATCONT.firstordercomponentnumber,[]),[2 1 3]),OCMATCONT.firstordercomponentnumber,[]);
            X0=reshape(permute(reshape(OCMATCONT.ZeroOrderCollocationValues(:,1:OCMATCONT.CollocationNumber+1).'*arc_coeff0,OCMATCONT.CollocationNumber+1,OCMATCONT.zeroordercomponentnumber,[]),[2 1 3]),OCMATCONT.zeroordercomponentnumber,[]);
            X1T=OCMATCONT.FirstOrderCollocationValues(:,OCMATCONT.CollocationNumber+2).'*arc_coeff1(:,end-OCMATCONT.firstordercomponentnumber+1:end);
            X0T=OCMATCONT.ZeroOrderCollocationValues(:,OCMATCONT.CollocationNumber+2).'*arc_coeff0(:,end-OCMATCONT.zeroordercomponentnumber+1:end);
            X(OCMATCONT.firstordercoordinate,ctr2+(1:((OCMATCONT.meshNumber(ii)-1)*(OCMATCONT.CollocationNumber+1)+1)))=[X1 X1T(:)];
            X(OCMATCONT.zeroordercoordinate,ctr2+(1:((OCMATCONT.meshNumber(ii)-1)*(OCMATCONT.CollocationNumber+1)+1)))=[X0 X0T(:)];
            ctr2=ctr2+((OCMATCONT.meshNumber(ii)-1)*(OCMATCONT.CollocationNumber+1)+1);
        end
    case 'grid'
        %X10=reshape(OCMATCONT.FirstOrderCollocationValues(:,1).'*coeff1,OCMATCONT.firstordercomponentnumber,[]);
        %X00=reshape(OCMATCONT.ZeroOrderCollocationValues(:,1).'*coeff0,OCMATCONT.zeroordercomponentnumber,[]);

        X=zeros(OCMATCONT.componentnumber,sum(OCMATCONT.meshNumber));
        ctr1=0;
        ctr0=0;
        ctr2=0;
        for ii=1:OCMATCONT.arcnumber
            arc_coeff1=coeff1(:,ctr1+(1:((OCMATCONT.meshNumber(ii)-1)*OCMATCONT.firstordercomponentnumber)));
            arc_coeff0=coeff0(:,ctr0+(1:((OCMATCONT.meshNumber(ii)-1)*OCMATCONT.zeroordercomponentnumber)));
            ctr1=ctr1+(OCMATCONT.meshNumber(ii)-1)*OCMATCONT.firstordercomponentnumber;
            ctr0=ctr0+(OCMATCONT.meshNumber(ii)-1)*OCMATCONT.zeroordercomponentnumber;
            X1=reshape(OCMATCONT.FirstOrderCollocationValues(:,1).'*arc_coeff1,OCMATCONT.firstordercomponentnumber,[]);
            X0=reshape(OCMATCONT.ZeroOrderCollocationValues(:,1).'*arc_coeff0,OCMATCONT.zeroordercomponentnumber,[]);
            X1T=OCMATCONT.FirstOrderCollocationValues(:,OCMATCONT.CollocationNumber+2).'*arc_coeff1(:,end-OCMATCONT.firstordercomponentnumber+1:end);
            X0T=OCMATCONT.ZeroOrderCollocationValues(:,OCMATCONT.CollocationNumber+2).'*arc_coeff0(:,end-OCMATCONT.zeroordercomponentnumber+1:end);
            X(OCMATCONT.firstordercoordinate,ctr2+(1:OCMATCONT.meshNumber(ii)))=[X1 X1T(:)];
            X(OCMATCONT.zeroordercoordinate,ctr2+(1:OCMATCONT.meshNumber(ii)))=[X0 X0T(:)];
            ctr2=ctr2+OCMATCONT.meshNumber(ii);
        end
end