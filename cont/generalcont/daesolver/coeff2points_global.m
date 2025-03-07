function X=coeff2points_global(tcolmesh,coeff,flag,varargin)
global OCMATCONT
switch flag
    case 'collocationgrid'
        % the components are computed als the polynomial values evaluated
        % at 0 rho1 ... rho_m
        coeff1=coeff(OCMATCONT.MESHDATA.firstorderidx);
        coeff0=coeff(OCMATCONT.MESHDATA.zeroorderidx);
        X=zeros(OCMATCONT.componentnumber,sum(OCMATCONT.MESHDATA.meshNumber-1)*(OCMATCONT.CollocationNumber+1)+OCMATCONT.arcnumber);
        ctr1=0;
        ctr0=0;
        ctr2=0;
        tmesh=tcolmesh(OCMATCONT.MESHDATA.tmeshidx);
        for ii=1:OCMATCONT.arcnumber
            H=repmat(diff(tmesh(OCMATCONT.MESHDATA.arcposition(1,ii):OCMATCONT.MESHDATA.arcposition(2,ii))),OCMATCONT.firstordercomponentnumber,1);   
            H=repmat(H(:).',OCMATCONT.CollocationNumber,1);
            arc_coeff1=coeff1(:,ctr1+(1:((OCMATCONT.MESHDATA.meshNumber(ii)-1)*OCMATCONT.firstordercomponentnumber)));
            arc_coeff1(2:OCMATCONT.CollocationNumber+1,:)=arc_coeff1(2:OCMATCONT.CollocationNumber+1,:).*H;
            arc_coeff0=coeff0(:,ctr0+(1:((OCMATCONT.MESHDATA.meshNumber(ii)-1)*OCMATCONT.zeroordercomponentnumber)));
            ctr1=ctr1+(OCMATCONT.MESHDATA.meshNumber(ii)-1)*OCMATCONT.firstordercomponentnumber;
            ctr0=ctr0+(OCMATCONT.MESHDATA.meshNumber(ii)-1)*OCMATCONT.zeroordercomponentnumber;
            X1=reshape(permute(reshape(OCMATCONT.FirstOrderCollocationValues(:,1:OCMATCONT.CollocationNumber+1).'*arc_coeff1,OCMATCONT.CollocationNumber+1,OCMATCONT.firstordercomponentnumber,[]),[2 1 3]),OCMATCONT.firstordercomponentnumber,[]);
            X0=reshape(permute(reshape(OCMATCONT.ZeroOrderCollocationValues(:,1:OCMATCONT.CollocationNumber+1).'*arc_coeff0,OCMATCONT.CollocationNumber+1,OCMATCONT.zeroordercomponentnumber,[]),[2 1 3]),OCMATCONT.zeroordercomponentnumber,[]);
            X1T=OCMATCONT.FirstOrderCollocationValues(:,OCMATCONT.CollocationNumber+2).'*arc_coeff1(:,end-OCMATCONT.firstordercomponentnumber+1:end);
            X0T=OCMATCONT.ZeroOrderCollocationValues(:,OCMATCONT.CollocationNumber+2).'*arc_coeff0(:,end-OCMATCONT.zeroordercomponentnumber+1:end);
            X(OCMATCONT.firstordercoordinate,ctr2+(1:((OCMATCONT.MESHDATA.meshNumber(ii)-1)*(OCMATCONT.CollocationNumber+1)+1)))=[X1 X1T(:)];
            X(OCMATCONT.zeroordercoordinate,ctr2+(1:((OCMATCONT.MESHDATA.meshNumber(ii)-1)*(OCMATCONT.CollocationNumber+1)+1)))=[X0 X0T(:)];
            ctr2=ctr2+((OCMATCONT.MESHDATA.meshNumber(ii)-1)*(OCMATCONT.CollocationNumber+1)+1);
        end
    case 'collocationgridmT'
        % the components are computed als the polynomial values evaluated
        % at 0 rho1 ... rho_m
        coeff1=coeff(OCMATCONT.MESHDATA.firstorderidx);
        coeff0=coeff(OCMATCONT.MESHDATA.zeroorderidx);
        X=zeros(OCMATCONT.componentnumber,sum(OCMATCONT.MESHDATA.meshNumber-1)*(OCMATCONT.CollocationNumber+1));
        ctr1=0;
        ctr0=0;
        ctr2=0;
        for ii=1:OCMATCONT.arcnumber
            H=repmat(diff(tcolmesh(OCMATCONT.MESHDATA.arcposition(1,ii):OCMATCONT.MESHDATA.arcposition(2,ii))),OCMATCONT.firstordercomponentnumber,1);   
            H=repmat(H(:).',OCMATCONT.CollocationNumber,1);
            arc_coeff1=coeff1(:,ctr1+(1:((OCMATCONT.MESHDATA.meshNumber(ii)-1)*OCMATCONT.firstordercomponentnumber)));
            arc_coeff1(2:OCMATCONT.CollocationNumber+1,:)=arc_coeff1(2:OCMATCONT.CollocationNumber+1,:).*H;
            arc_coeff0=coeff0(:,ctr0+(1:((OCMATCONT.MESHDATA.meshNumber(ii)-1)*OCMATCONT.zeroordercomponentnumber)));
            ctr1=ctr1+(OCMATCONT.MESHDATA.meshNumber(ii)-1)*OCMATCONT.firstordercomponentnumber;
            ctr0=ctr0+(OCMATCONT.MESHDATA.meshNumber(ii)-1)*OCMATCONT.zeroordercomponentnumber;
            X(OCMATCONT.firstordercoordinate,ctr2+(1:((OCMATCONT.MESHDATA.meshNumber(ii)-1)*(OCMATCONT.CollocationNumber+1))))=reshape(permute(reshape(OCMATCONT.FirstOrderCollocationValues(:,1:OCMATCONT.CollocationNumber+1).'*arc_coeff1,OCMATCONT.CollocationNumber+1,OCMATCONT.firstordercomponentnumber,[]),[2 1 3]),OCMATCONT.firstordercomponentnumber,[]);
            X(OCMATCONT.zeroordercoordinate,ctr2+(1:((OCMATCONT.MESHDATA.meshNumber(ii)-1)*(OCMATCONT.CollocationNumber+1))))=reshape(permute(reshape(OCMATCONT.ZeroOrderCollocationValues(:,1:OCMATCONT.CollocationNumber+1).'*arc_coeff0,OCMATCONT.CollocationNumber+1,OCMATCONT.zeroordercomponentnumber,[]),[2 1 3]),OCMATCONT.zeroordercomponentnumber,[]);
            ctr2=ctr2+((OCMATCONT.MESHDATA.meshNumber(ii)-1)*(OCMATCONT.CollocationNumber+1)+1);
        end
    case 'collocationgrid_'
        % the components are computed als the polynomial values evaluated
        % at rho1 ... rho_m 1

        coeff1=coeff(OCMATCONT.MESHDATA.firstorderidx);
        coeff0=coeff(OCMATCONT.MESHDATA.zeroorderidx);
        
        X=zeros(OCMATCONT.componentnumber,sum(OCMATCONT.MESHDATA.meshNumber-1)*(OCMATCONT.CollocationNumber+1)+OCMATCONT.arcnumber);
        ctr1=0;
        ctr0=0;
        ctr2=0;
        for ii=1:OCMATCONT.arcnumber
            H=repmat(diff(tcolmesh(OCMATCONT.MESHDATA.arcposition(1,ii):OCMATCONT.MESHDATA.arcposition(2,ii))),OCMATCONT.firstordercomponentnumber,1);   
            H=repmat(H(:).',OCMATCONT.CollocationNumber,1);
            arc_coeff1=coeff1(:,ctr1+(1:((OCMATCONT.MESHDATA.meshNumber(ii)-1)*OCMATCONT.firstordercomponentnumber)));
            arc_coeff1(2:OCMATCONT.CollocationNumber+1,:)=arc_coeff1(2:OCMATCONT.CollocationNumber+1,:).*H;
            arc_coeff0=coeff0(:,ctr0+(1:((OCMATCONT.MESHDATA.meshNumber(ii)-1)*OCMATCONT.zeroordercomponentnumber)));
            ctr1=ctr1+(OCMATCONT.MESHDATA.meshNumber(ii)-1)*OCMATCONT.firstordercomponentnumber;
            ctr0=ctr0+(OCMATCONT.MESHDATA.meshNumber(ii)-1)*OCMATCONT.zeroordercomponentnumber;
            X1=reshape(permute(reshape(OCMATCONT.FirstOrderCollocationValues(:,2:OCMATCONT.CollocationNumber+2).'*arc_coeff1,OCMATCONT.CollocationNumber+1,OCMATCONT.firstordercomponentnumber,[]),[2 1 3]),OCMATCONT.firstordercomponentnumber,[]);
            X0=reshape(permute(reshape(OCMATCONT.ZeroOrderCollocationValues(:,2:OCMATCONT.CollocationNumber+2).'*arc_coeff0,OCMATCONT.CollocationNumber+1,OCMATCONT.zeroordercomponentnumber,[]),[2 1 3]),OCMATCONT.zeroordercomponentnumber,[]);
            X10=OCMATCONT.FirstOrderCollocationValues(:,1).'*arc_coeff1(:,1:OCMATCONT.firstordercomponentnumber);
            X00=OCMATCONT.ZeroOrderCollocationValues(:,1).'*arc_coeff0(:,1:OCMATCONT.zeroordercomponentnumber);
            X(OCMATCONT.firstordercoordinate,ctr2+(1:((OCMATCONT.MESHDATA.meshNumber(ii)-1)*(OCMATCONT.CollocationNumber+1)+1)))=[X10(:) X1];
            X(OCMATCONT.zeroordercoordinate,ctr2+(1:((OCMATCONT.MESHDATA.meshNumber(ii)-1)*(OCMATCONT.CollocationNumber+1)+1)))=[X00(:) X0];
            ctr2=ctr2+((OCMATCONT.MESHDATA.meshNumber(ii)-1)*(OCMATCONT.CollocationNumber+1)+1);
        end
    case 'grid'

        tmesh=tcolmesh(OCMATCONT.MESHDATA.tmeshidx);
        coeff1=coeff(OCMATCONT.MESHDATA.firstorderidx);
        coeff0=coeff(OCMATCONT.MESHDATA.zeroorderidx);

        X=zeros(OCMATCONT.componentnumber,sum(OCMATCONT.MESHDATA.meshNumber));
        ctr1=0;
        ctr0=0;
        ctr2=0;
        for ii=1:OCMATCONT.arcnumber
            H=repmat(diff(tmesh(OCMATCONT.MESHDATA.arcposition(1,ii):OCMATCONT.MESHDATA.arcposition(2,ii))),OCMATCONT.firstordercomponentnumber,1);   
            H=repmat(H(:).',OCMATCONT.CollocationNumber,1);
            arc_coeff1=coeff1(:,ctr1+(1:((OCMATCONT.MESHDATA.meshNumber(ii)-1)*OCMATCONT.firstordercomponentnumber)));
            arc_coeff1(2:OCMATCONT.CollocationNumber+1,:)=arc_coeff1(2:OCMATCONT.CollocationNumber+1,:).*H;
            arc_coeff0=coeff0(:,ctr0+(1:((OCMATCONT.MESHDATA.meshNumber(ii)-1)*OCMATCONT.zeroordercomponentnumber)));
            ctr1=ctr1+(OCMATCONT.MESHDATA.meshNumber(ii)-1)*OCMATCONT.firstordercomponentnumber;
            ctr0=ctr0+(OCMATCONT.MESHDATA.meshNumber(ii)-1)*OCMATCONT.zeroordercomponentnumber;
            X1=reshape(OCMATCONT.FirstOrderCollocationValues(:,1).'*arc_coeff1,OCMATCONT.firstordercomponentnumber,[]);
            X0=reshape(OCMATCONT.ZeroOrderCollocationValues(:,1).'*arc_coeff0,OCMATCONT.zeroordercomponentnumber,[]);
            X1T=OCMATCONT.FirstOrderCollocationValues(:,OCMATCONT.CollocationNumber+2).'*arc_coeff1(:,end-OCMATCONT.firstordercomponentnumber+1:end);
            X0T=OCMATCONT.ZeroOrderCollocationValues(:,OCMATCONT.CollocationNumber+2).'*arc_coeff0(:,end-OCMATCONT.zeroordercomponentnumber+1:end);
            X(OCMATCONT.firstordercoordinate,ctr2+(1:OCMATCONT.MESHDATA.meshNumber(ii)))=[X1 X1T(:)];
            X(OCMATCONT.zeroordercoordinate,ctr2+(1:OCMATCONT.MESHDATA.meshNumber(ii)))=[X0 X0T(:)];
            ctr2=ctr2+OCMATCONT.MESHDATA.meshNumber(ii);
        end
    case 'collocationgrid1_2'
        % the components are computed als the polynomial values evaluated
        % at 0 rho1 ... rho_m
        tmesh=tcolmesh(OCMATCONT.MESHDATA1_2.tmeshidx);
        coeff1=coeff(OCMATCONT.MESHDATA1_2.firstorderidx);
        coeff0=coeff(OCMATCONT.MESHDATA1_2.zeroorderidx);

        X=zeros(OCMATCONT.componentnumber,sum(OCMATCONT.MESHDATA1_2.meshNumber-1)*(OCMATCONT.CollocationNumber+1)+OCMATCONT.arcnumber);
        ctr1=0;
        ctr0=0;
        ctr2=0;
        for ii=1:OCMATCONT.arcnumber
            H=repmat(diff(tmesh(OCMATCONT.MESHDATA1_2.arcposition(1,ii):OCMATCONT.MESHDATA1_2.arcposition(2,ii))),OCMATCONT.firstordercomponentnumber,1);
            H=repmat(H(:).',OCMATCONT.CollocationNumber,1);
            arc_coeff1=coeff1(:,ctr1+(1:((OCMATCONT.MESHDATA1_2.meshNumber(ii)-1)*OCMATCONT.firstordercomponentnumber)));
            arc_coeff1(2:OCMATCONT.CollocationNumber+1,:)=arc_coeff1(2:OCMATCONT.CollocationNumber+1,:).*H;
            arc_coeff0=coeff0(:,ctr0+(1:((OCMATCONT.MESHDATA1_2.meshNumber(ii)-1)*OCMATCONT.zeroordercomponentnumber)));
            ctr1=ctr1+(OCMATCONT.MESHDATA1_2.meshNumber(ii)-1)*OCMATCONT.firstordercomponentnumber;
            ctr0=ctr0+(OCMATCONT.MESHDATA1_2.meshNumber(ii)-1)*OCMATCONT.zeroordercomponentnumber;

            X1=reshape(permute(reshape(OCMATCONT.FirstOrderCollocationValues(:,1:OCMATCONT.CollocationNumber+1).'*arc_coeff1,OCMATCONT.CollocationNumber+1,OCMATCONT.firstordercomponentnumber,[]),[2 1 3]),OCMATCONT.firstordercomponentnumber,[]);
            X0=reshape(permute(reshape(OCMATCONT.ZeroOrderCollocationValues(:,1:OCMATCONT.CollocationNumber+1).'*arc_coeff0,OCMATCONT.CollocationNumber+1,OCMATCONT.zeroordercomponentnumber,[]),[2 1 3]),OCMATCONT.zeroordercomponentnumber,[]);
            X1T=OCMATCONT.FirstOrderCollocationValues(:,OCMATCONT.CollocationNumber+2).'*arc_coeff1(:,end-OCMATCONT.firstordercomponentnumber+1:end);
            X0T=OCMATCONT.ZeroOrderCollocationValues(:,OCMATCONT.CollocationNumber+2).'*arc_coeff0(:,end-OCMATCONT.zeroordercomponentnumber+1:end);
            X(OCMATCONT.firstordercoordinate,ctr2+(1:((OCMATCONT.MESHDATA1_2.meshNumber(ii)-1)*(OCMATCONT.CollocationNumber+1)+1)))=[X1 X1T(:)];
            X(OCMATCONT.zeroordercoordinate,ctr2+(1:((OCMATCONT.MESHDATA1_2.meshNumber(ii)-1)*(OCMATCONT.CollocationNumber+1)+1)))=[X0 X0T(:)];
            ctr2=ctr2+((OCMATCONT.MESHDATA1_2.meshNumber(ii)-1)*(OCMATCONT.CollocationNumber+1)+1);
        end
    case 'collocationgrid1_2_'
        % the components are computed als the polynomial values evaluated
        % at rho1 ... rho_m 1

        coeff1=coeff(OCMATCONT.MESHDATA1_2.firstorderidx);
        coeff0=coeff(OCMATCONT.MESHDATA1_2.zeroorderidx);
        
        X=zeros(OCMATCONT.componentnumber,sum(OCMATCONT.MESHDATA1_2.meshNumber-1)*(OCMATCONT.CollocationNumber+1)+OCMATCONT.arcnumber);
        ctr1=0;
        ctr0=0;
        ctr2=0;
        for ii=1:OCMATCONT.arcnumber
            H=repmat(diff(tcolmesh(OCMATCONT.MESHDATA1_2.arcposition(1,ii):OCMATCONT.MESHDATA1_2.arcposition(2,ii))),OCMATCONT.firstordercomponentnumber,1);   
            H=repmat(H(:).',OCMATCONT.CollocationNumber,1);
            arc_coeff1=coeff1(:,ctr1+(1:((OCMATCONT.MESHDATA1_2.meshNumber(ii)-1)*OCMATCONT.firstordercomponentnumber)));
            arc_coeff1(2:OCMATCONT.CollocationNumber+1,:)=arc_coeff1(2:OCMATCONT.CollocationNumber+1,:).*H;
            arc_coeff0=coeff0(:,ctr0+(1:((OCMATCONT.MESHDATA1_2.meshNumber(ii)-1)*OCMATCONT.zeroordercomponentnumber)));
            ctr1=ctr1+(OCMATCONT.MESHDATA1_2.meshNumber(ii)-1)*OCMATCONT.firstordercomponentnumber;
            ctr0=ctr0+(OCMATCONT.MESHDATA1_2.meshNumber(ii)-1)*OCMATCONT.zeroordercomponentnumber;
            X1=reshape(permute(reshape(OCMATCONT.FirstOrderCollocationValues(:,2:OCMATCONT.CollocationNumber+2).'*arc_coeff1,OCMATCONT.CollocationNumber+1,OCMATCONT.firstordercomponentnumber,[]),[2 1 3]),OCMATCONT.firstordercomponentnumber,[]);
            X0=reshape(permute(reshape(OCMATCONT.ZeroOrderCollocationValues(:,2:OCMATCONT.CollocationNumber+2).'*arc_coeff0,OCMATCONT.CollocationNumber+1,OCMATCONT.zeroordercomponentnumber,[]),[2 1 3]),OCMATCONT.zeroordercomponentnumber,[]);
            X10=OCMATCONT.FirstOrderCollocationValues(:,1).'*arc_coeff1(:,1:OCMATCONT.firstordercomponentnumber);
            X00=OCMATCONT.ZeroOrderCollocationValues(:,1).'*arc_coeff0(:,1:OCMATCONT.zeroordercomponentnumber);
            X(OCMATCONT.firstordercoordinate,ctr2+(1:((OCMATCONT.MESHDATA1_2.meshNumber(ii)-1)*(OCMATCONT.CollocationNumber+1)+1)))=[X10(:) X1];
            X(OCMATCONT.zeroordercoordinate,ctr2+(1:((OCMATCONT.MESHDATA1_2.meshNumber(ii)-1)*(OCMATCONT.CollocationNumber+1)+1)))=[X00(:) X0];
            ctr2=ctr2+((OCMATCONT.MESHDATA1_2.meshNumber(ii)-1)*(OCMATCONT.CollocationNumber+1)+1);
        end
    case 'grid1_2'

        tmesh=tcolmesh(OCMATCONT.MESHDATA1_2.tmeshidx);
        coeff1=coeff(OCMATCONT.MESHDATA1_2.firstorderidx);
        coeff0=coeff(OCMATCONT.MESHDATA1_2.zeroorderidx);

        X=zeros(OCMATCONT.componentnumber,sum(OCMATCONT.MESHDATA1_2.meshNumber));
        ctr1=0;
        ctr0=0;
        ctr2=0;
        for ii=1:OCMATCONT.arcnumber
            H=repmat(diff(tmesh(OCMATCONT.MESHDATA1_2.arcposition(1,ii):OCMATCONT.MESHDATA1_2.arcposition(2,ii))),OCMATCONT.firstordercomponentnumber,1);   
            H=repmat(H(:).',OCMATCONT.CollocationNumber,1);
            arc_coeff1=coeff1(:,ctr1+(1:((OCMATCONT.MESHDATA1_2.meshNumber(ii)-1)*OCMATCONT.firstordercomponentnumber)));
            arc_coeff1(2:OCMATCONT.CollocationNumber+1,:)=arc_coeff1(2:OCMATCONT.CollocationNumber+1,:).*H;
            arc_coeff0=coeff0(:,ctr0+(1:((OCMATCONT.MESHDATA1_2.meshNumber(ii)-1)*OCMATCONT.zeroordercomponentnumber)));
            ctr1=ctr1+(OCMATCONT.MESHDATA1_2.meshNumber(ii)-1)*OCMATCONT.firstordercomponentnumber;
            ctr0=ctr0+(OCMATCONT.MESHDATA1_2.meshNumber(ii)-1)*OCMATCONT.zeroordercomponentnumber;
            X1=reshape(OCMATCONT.FirstOrderCollocationValues(:,1).'*arc_coeff1,OCMATCONT.firstordercomponentnumber,[]);
            X0=reshape(OCMATCONT.ZeroOrderCollocationValues(:,1).'*arc_coeff0,OCMATCONT.zeroordercomponentnumber,[]);
            X1T=OCMATCONT.FirstOrderCollocationValues(:,OCMATCONT.CollocationNumber+2).'*arc_coeff1(:,end-OCMATCONT.firstordercomponentnumber+1:end);
            X0T=OCMATCONT.ZeroOrderCollocationValues(:,OCMATCONT.CollocationNumber+2).'*arc_coeff0(:,end-OCMATCONT.zeroordercomponentnumber+1:end);
            X(OCMATCONT.firstordercoordinate,ctr2+(1:OCMATCONT.MESHDATA1_2.meshNumber(ii)))=[X1 X1T(:)];
            X(OCMATCONT.zeroordercoordinate,ctr2+(1:OCMATCONT.MESHDATA1_2.meshNumber(ii)))=[X0 X0T(:)];
            ctr2=ctr2+OCMATCONT.MESHDATA1_2.meshNumber(ii);
        end
        
    case 'collocationgrid2collocationgrid1_2'
        % evaluates coeff on coarse grid at fine grid (fine grid = coarse
        % grid + midpoints and collocation points

        coeff1=coeff(OCMATCONT.MESHDATA.firstorderidx);
        coeff0=coeff(OCMATCONT.MESHDATA.zeroorderidx);
        X=zeros(OCMATCONT.componentnumber,sum(OCMATCONT.MESHDATA1_2.meshNumber-1)*(OCMATCONT.CollocationNumber+1)+OCMATCONT.arcnumber);
        ctr1=0;
        ctr0=0;
        ctr2=0;
        for ii=1:OCMATCONT.arcnumber
            H=repmat(diff(tcolmesh(OCMATCONT.MESHDATA.arcposition(1,ii):OCMATCONT.MESHDATA.arcposition(2,ii))),OCMATCONT.firstordercomponentnumber,1);   
            H=repmat(H(:).',OCMATCONT.CollocationNumber,1);
            arc_coeff1=coeff1(:,ctr1+(1:((OCMATCONT.MESHDATA.meshNumber(ii)-1)*OCMATCONT.firstordercomponentnumber)));
            arc_coeff1(2:OCMATCONT.CollocationNumber+1,:)=arc_coeff1(2:OCMATCONT.CollocationNumber+1,:).*H;
            arc_coeff0=coeff0(:,ctr0+(1:((OCMATCONT.MESHDATA.meshNumber(ii)-1)*OCMATCONT.zeroordercomponentnumber)));
            ctr1=ctr1+(OCMATCONT.MESHDATA.meshNumber(ii)-1)*OCMATCONT.firstordercomponentnumber;
            ctr0=ctr0+(OCMATCONT.MESHDATA.meshNumber(ii)-1)*OCMATCONT.zeroordercomponentnumber;
            X1=reshape(permute(reshape(OCMATCONT.FirstOrderCollocationValues1_2(:,1:2*(OCMATCONT.CollocationNumber+1)).'*arc_coeff1,2*(OCMATCONT.CollocationNumber+1),OCMATCONT.firstordercomponentnumber,[]),[2 1 3]),OCMATCONT.firstordercomponentnumber,[]);
            X0=reshape(permute(reshape(OCMATCONT.ZeroOrderCollocationValues1_2(:,1:2*(OCMATCONT.CollocationNumber+1)).'*arc_coeff0,2*(OCMATCONT.CollocationNumber+1),OCMATCONT.zeroordercomponentnumber,[]),[2 1 3]),OCMATCONT.zeroordercomponentnumber,[]);
            X1T=OCMATCONT.FirstOrderCollocationValues(:,OCMATCONT.CollocationNumber+2).'*arc_coeff1(:,end-OCMATCONT.firstordercomponentnumber+1:end);
            X0T=OCMATCONT.ZeroOrderCollocationValues(:,OCMATCONT.CollocationNumber+2).'*arc_coeff0(:,end-OCMATCONT.zeroordercomponentnumber+1:end);
            X(OCMATCONT.firstordercoordinate,ctr2+(1:((OCMATCONT.MESHDATA.meshNumber(ii)-1)*(2*(OCMATCONT.CollocationNumber+1))+1)))=[X1 X1T(:) ];
            X(OCMATCONT.zeroordercoordinate,ctr2+(1:((OCMATCONT.MESHDATA.meshNumber(ii)-1)*(2*(OCMATCONT.CollocationNumber+1))+1)))=[X0 X0T(:)];
            ctr2=ctr2+((OCMATCONT.MESHDATA.meshNumber(ii)-1)*(2*(OCMATCONT.CollocationNumber+1))+1);
        end
        
    case 'collocationgrid2collocationgrid1_2_'
        % evaluates coeff on coarse grid at fine grid (fine grid = coarse
        % grid + midpoints and collocation points

        coeff1=coeff(OCMATCONT.MESHDATA.firstorderidx);
        coeff0=coeff(OCMATCONT.MESHDATA.zeroorderidx);
        X=zeros(OCMATCONT.componentnumber,sum(OCMATCONT.MESHDATA1_2.meshNumber-1)*(OCMATCONT.CollocationNumber+1)+OCMATCONT.arcnumber);
        ctr1=0;
        ctr0=0;
        ctr2=0;
        for ii=1:OCMATCONT.arcnumber
            H=repmat(diff(tcolmesh(OCMATCONT.MESHDATA.arcposition(1,ii):OCMATCONT.MESHDATA.arcposition(2,ii))),OCMATCONT.firstordercomponentnumber,1);   
            H=repmat(H(:).',OCMATCONT.CollocationNumber,1);
            arc_coeff1=coeff1(:,ctr1+(1:((OCMATCONT.MESHDATA.meshNumber(ii)-1)*OCMATCONT.firstordercomponentnumber)));
            arc_coeff1(2:OCMATCONT.CollocationNumber+1,:)=arc_coeff1(2:OCMATCONT.CollocationNumber+1,:).*H;
            arc_coeff0=coeff0(:,ctr0+(1:((OCMATCONT.MESHDATA.meshNumber(ii)-1)*OCMATCONT.zeroordercomponentnumber)));
            ctr1=ctr1+(OCMATCONT.MESHDATA.meshNumber(ii)-1)*OCMATCONT.firstordercomponentnumber;
            ctr0=ctr0+(OCMATCONT.MESHDATA.meshNumber(ii)-1)*OCMATCONT.zeroordercomponentnumber;
            X1=reshape(permute(reshape(OCMATCONT.FirstOrderCollocationValues1_2(:,2:2*(OCMATCONT.CollocationNumber+1)+1).'*arc_coeff1,2*(OCMATCONT.CollocationNumber+1),OCMATCONT.firstordercomponentnumber,[]),[2 1 3]),OCMATCONT.firstordercomponentnumber,[]);
            X0=reshape(permute(reshape(OCMATCONT.ZeroOrderCollocationValues1_2(:,2:2*(OCMATCONT.CollocationNumber+1)+1).'*arc_coeff0,2*(OCMATCONT.CollocationNumber+1),OCMATCONT.zeroordercomponentnumber,[]),[2 1 3]),OCMATCONT.zeroordercomponentnumber,[]);
            X10=OCMATCONT.FirstOrderCollocationValues1_2(:,1).'*arc_coeff1(:,1:OCMATCONT.firstordercomponentnumber);
            X00=OCMATCONT.ZeroOrderCollocationValues1_2(:,1).'*arc_coeff0(:,1:OCMATCONT.zeroordercomponentnumber);
            X(OCMATCONT.firstordercoordinate,ctr2+(1:((OCMATCONT.MESHDATA.meshNumber(ii)-1)*(2*(OCMATCONT.CollocationNumber+1))+1)))=[X10(:) X1];
            X(OCMATCONT.zeroordercoordinate,ctr2+(1:((OCMATCONT.MESHDATA.meshNumber(ii)-1)*(2*(OCMATCONT.CollocationNumber+1))+1)))=[X00(:) X0];
            ctr2=ctr2+((OCMATCONT.MESHDATA.meshNumber(ii)-1)*(2*(OCMATCONT.CollocationNumber+1))+1);
        end
        
    case 'collocationgrid1_22collocationgrid'
        % evaluates coeff on fine grid at coarse grid (fine grid = coarse
        % grid + midpoints and collocation points
        tmesh=tcolmesh(OCMATCONT.MESHDATA1_2.tmeshidx);
        coeff1=coeff(OCMATCONT.MESHDATA1_2.firstorderidx);
        coeff0=coeff(OCMATCONT.MESHDATA1_2.zeroorderidx);
        X=zeros(OCMATCONT.componentnumber,sum(OCMATCONT.MESHDATA.meshNumber-1)*(OCMATCONT.NumberPoints2_1-1)+OCMATCONT.arcnumber);
        ctr1=0;
        ctr0=0;
        ctr2=0;
        for ii=1:OCMATCONT.arcnumber
            H=repmat(diff(tmesh(OCMATCONT.MESHDATA1_2.arcposition(1,ii):OCMATCONT.MESHDATA1_2.arcposition(2,ii))),OCMATCONT.firstordercomponentnumber,1);   
            H=repmat(H(:).',OCMATCONT.CollocationNumber,1);
%             keepidx=2*(0:((OCMATCONT.MESHDATA1_2.meshNumber(ii)-2)/2))*(OCMATCONT.NumberPoints2_1-1);
%             keepidx=[repmat(keepidx(:),1,numel(OCMATCONT.IndexLeft2_1))+repmat(OCMATCONT.IndexLeft2_1,length(keepidx),1) repmat(keepidx(:),1,numel(OCMATCONT.IndexRight2_1)-1)+OCMATCONT.NumberPoints2_1-1+repmat(OCMATCONT.IndexRight2_1(1:end-1),length(keepidx),1)];
%             keepidx=keepidx.';
%             keepidx=[keepidx(:).' (OCMATCONT.MESHDATA1_2.meshNumber(ii)-1)*(OCMATCONT.NumberPoints2_1-1)+1];
            arc_coeff1=coeff1(:,ctr1+(1:((OCMATCONT.MESHDATA1_2.meshNumber(ii)-1)*OCMATCONT.firstordercomponentnumber)));
            arc_coeff1(2:OCMATCONT.CollocationNumber+1,:)=arc_coeff1(2:OCMATCONT.CollocationNumber+1,:).*H;
            arc_coeff0=coeff0(:,ctr0+(1:((OCMATCONT.MESHDATA1_2.meshNumber(ii)-1)*OCMATCONT.zeroordercomponentnumber)));
            ctr1=ctr1+(OCMATCONT.MESHDATA1_2.meshNumber(ii)-1)*OCMATCONT.firstordercomponentnumber;
            ctr0=ctr0+(OCMATCONT.MESHDATA1_2.meshNumber(ii)-1)*OCMATCONT.zeroordercomponentnumber;
            X1=reshape(permute(reshape(OCMATCONT.FirstOrderCollocationValues2_1(:,1:OCMATCONT.NumberPoints2_1-1).'*arc_coeff1,OCMATCONT.NumberPoints2_1-1,OCMATCONT.firstordercomponentnumber,[]),[2 1 3]),OCMATCONT.firstordercomponentnumber,[]);
            X0=reshape(permute(reshape(OCMATCONT.ZeroOrderCollocationValues2_1(:,1:OCMATCONT.NumberPoints2_1-1).'*arc_coeff0,OCMATCONT.NumberPoints2_1-1,OCMATCONT.zeroordercomponentnumber,[]),[2 1 3]),OCMATCONT.zeroordercomponentnumber,[]);
            X1T=OCMATCONT.FirstOrderCollocationValues2_1(:,end).'*arc_coeff1(:,end-OCMATCONT.firstordercomponentnumber+1:end);
            X0T=OCMATCONT.ZeroOrderCollocationValues2_1(:,end).'*arc_coeff0(:,end-OCMATCONT.zeroordercomponentnumber+1:end);
            Xtmp=[[X1 X1T(:)];[X0 X0T(:)]];
            X(1:OCMATCONT.componentnumber,ctr2+(1:((OCMATCONT.MESHDATA.meshNumber(ii)-1)*(OCMATCONT.NumberPoints2_1-1)+1)))=Xtmp([OCMATCONT.firstordercoordinate;OCMATCONT.zeroordercoordinate],OCMATCONT.MESHDATA1_2.Matrix(ii).KeepIdx2_1);
            ctr2=ctr2+((OCMATCONT.MESHDATA.meshNumber(ii)-1)*(OCMATCONT.NumberPoints2_1-1)+1);
        end
        
    case 'collocationgrid1_22collocationgrid_'
        % evaluates coeff on fine grid at coarse grid (fine grid = coarse
        % grid + midpoints and collocation points
        tmesh=tcolmesh(OCMATCONT.MESHDATA1_2.tmeshidx);
        coeff1=coeff(OCMATCONT.MESHDATA1_2.firstorderidx);
        coeff0=coeff(OCMATCONT.MESHDATA1_2.zeroorderidx);
        X=zeros(OCMATCONT.componentnumber,sum(OCMATCONT.MESHDATA.meshNumber-1)*(OCMATCONT.NumberPoints2_1-1)+OCMATCONT.arcnumber);
        ctr1=0;
        ctr0=0;
        ctr2=0;
        for ii=1:OCMATCONT.arcnumber
            H=repmat(diff(tmesh(OCMATCONT.MESHDATA1_2.arcposition(1,ii):OCMATCONT.MESHDATA1_2.arcposition(2,ii))),OCMATCONT.firstordercomponentnumber,1);   
            H=repmat(H(:).',OCMATCONT.CollocationNumber,1);
%             keepidx=2*(0:((OCMATCONT.MESHDATA1_2.meshNumber(ii)-2)/2))*(OCMATCONT.NumberPoints2_1-1);
%             keepidx=[repmat(keepidx(:),1,numel(OCMATCONT.IndexLeft2_1))+repmat(OCMATCONT.IndexLeft2_1,length(keepidx),1) repmat(keepidx(:),1,numel(OCMATCONT.IndexRight2_1)-1)+OCMATCONT.NumberPoints2_1-1+repmat(OCMATCONT.IndexRight2_1(1:end-1),length(keepidx),1)];
%             keepidx=keepidx.';
%             keepidx=[keepidx(:).' (OCMATCONT.MESHDATA1_2.meshNumber(ii)-1)*(OCMATCONT.NumberPoints2_1-1)+1];
            arc_coeff1=coeff1(:,ctr1+(1:((OCMATCONT.MESHDATA1_2.meshNumber(ii)-1)*OCMATCONT.firstordercomponentnumber)));
            arc_coeff1(2:OCMATCONT.CollocationNumber+1,:)=arc_coeff1(2:OCMATCONT.CollocationNumber+1,:).*H;
            arc_coeff0=coeff0(:,ctr0+(1:((OCMATCONT.MESHDATA1_2.meshNumber(ii)-1)*OCMATCONT.zeroordercomponentnumber)));
            ctr1=ctr1+(OCMATCONT.MESHDATA1_2.meshNumber(ii)-1)*OCMATCONT.firstordercomponentnumber;
            ctr0=ctr0+(OCMATCONT.MESHDATA1_2.meshNumber(ii)-1)*OCMATCONT.zeroordercomponentnumber;
            X1=reshape(permute(reshape(OCMATCONT.FirstOrderCollocationValues2_1(:,2:OCMATCONT.NumberPoints2_1).'*arc_coeff1,OCMATCONT.NumberPoints2_1-1,OCMATCONT.firstordercomponentnumber,[]),[2 1 3]),OCMATCONT.firstordercomponentnumber,[]);
            X0=reshape(permute(reshape(OCMATCONT.ZeroOrderCollocationValues2_1(:,2:OCMATCONT.NumberPoints2_1).'*arc_coeff0,OCMATCONT.NumberPoints2_1-1,OCMATCONT.zeroordercomponentnumber,[]),[2 1 3]),OCMATCONT.zeroordercomponentnumber,[]);
            X10=OCMATCONT.FirstOrderCollocationValues2_1(:,1).'*arc_coeff1(:,1:OCMATCONT.firstordercomponentnumber);
            X00=OCMATCONT.ZeroOrderCollocationValues2_1(:,1).'*arc_coeff0(:,1:OCMATCONT.zeroordercomponentnumber);
            Xtmp=[[X10(:) X1];[X00(:) X0]];
            %X(1:OCMATCONT.componentnumber,ctr2+(1:((OCMATCONT.MESHDATA.meshNumber(ii)-1)*(OCMATCONT.NumberPoints2_1-1)+1)))=Xtmp(:,keepidx);
            X(1:OCMATCONT.componentnumber,ctr2+(1:((OCMATCONT.MESHDATA.meshNumber(ii)-1)*(OCMATCONT.NumberPoints2_1-1)+1)))=Xtmp([OCMATCONT.firstordercoordinate;OCMATCONT.zeroordercoordinate],OCMATCONT.MESHDATA1_2.Matrix(ii).KeepIdx2_1);
            ctr2=ctr2+((OCMATCONT.MESHDATA.meshNumber(ii)-1)*(OCMATCONT.NumberPoints2_1-1)+1);
        end
        
    case 'general'
        % coefficients based on data stored in MESHDATA
        if nargin==2
            return
        end
        t=varargin{1};
        tmesh=tcolmesh(OCMATCONT.MESHDATA.tmeshidx);
        X=zeros(OCMATCONT.componentnumber,length(t));
        coeff1=coeff(OCMATCONT.MESHDATA.firstorderidx);
        coeff0=coeff(OCMATCONT.MESHDATA.zeroorderidx);
        for ii=1:OCMATCONT.arcnumber
            h=diff(tmesh(OCMATCONT.MESHDATA.arcposition(1,ii):OCMATCONT.MESHDATA.arcposition(2,ii)));
            for jj=1:OCMATCONT.MESHDATA.meshNumber-1
                if jj<OCMATCONT.MESHDATA.meshNumber-1
                    % exclude the right side of the interval
                    tidx=find(tmesh(jj)<=t & t<tmesh(jj+1));
                else
                    % include the right side of the interval
                    tidx=find(tmesh(jj)<=t & t<=tmesh(jj+1));
                end
                normalizedt=(t(tidx)-tmesh(jj))/h(jj);
                coeff1part=coeff1(:,(jj-1)*OCMATCONT.firstordercomponentnumber+(1:OCMATCONT.firstordercomponentnumber));
                coeff0part=coeff0(:,(jj-1)*OCMATCONT.zeroordercomponentnumber+(1:OCMATCONT.zeroordercomponentnumber));
                psival1=ones(OCMATCONT.CollocationNumber+1,length(tidx));
                psival0=zeros(OCMATCONT.CollocationNumber,length(tidx));
                for m=1:OCMATCONT.CollocationNumber
                    psival1(m+1,:)=h(jj)*polyval(OCMATCONT.FirstOrderCollocationPolynomial(m,:),normalizedt);
                    psival0(m,:)=polyval(OCMATCONT.ZeroOrderCollocationPolynomial(m,:),normalizedt);
                end
                X(OCMATCONT.firstordercoordinate,tidx)=coeff1part.'*psival1;
                X(OCMATCONT.zeroordercoordinate,tidx)=coeff0part.'*psival0;
            end
        end
        
    case 'general1_2'
        % coefficients based on data stored in MESHDATA1_2 
        if nargin==2
            return
        end
        t=varargin{1};
        tmesh=tcolmesh(OCMATCONT.MESHDATA1_2.tmeshidx);
        coeff1=coeff(OCMATCONT.MESHDATA1_2.firstorderidx);
        coeff0=coeff(OCMATCONT.MESHDATA1_2.zeroorderidx);
        for ii=1:OCMATCONT.arcnumber
            h=diff(tmesh(OCMATCONT.MESHDATA1_2.arcposition(1,ii):OCMATCONT.MESHDATA1_2.arcposition(2,ii)));
            X=zeros(OCMATCONT.componentnumber,length(tmesh));
            for jj=1:OCMATCONT.MESHDATA1_2.meshNumber-1
                if jj<OCMATCONT.MESHDATA1_2.meshNumber-1
                    % exclude the right side of the interval
                    tidx=find(tmesh(jj)<=t & t<tmesh(jj+1));
                else
                    % include the right side of the interval
                    tidx=find(tmesh(jj)<=t & t<=tmesh(jj+1));
                end
                normalizedt=(t(tidx)-tmesh(jj))/h(jj);
                coeff1part=coeff1(:,(jj-1)*OCMATCONT.firstordercomponentnumber+(1:OCMATCONT.firstordercomponentnumber));
                coeff0part=coeff0(:,(jj-1)*OCMATCONT.zeroordercomponentnumber+(1:OCMATCONT.zeroordercomponentnumber));
                psival1=ones(OCMATCONT.CollocationNumber+1,length(tidx));
                psival0=zeros(OCMATCONT.CollocationNumber,length(tidx));
                for m=1:OCMATCONT.CollocationNumber
                    psival1(m+1,:)=h(jj)*polyval(OCMATCONT.FirstOrderCollocationPolynomial(m,:),normalizedt);
                    psival0(m,:)=polyval(OCMATCONT.ZeroOrderCollocationPolynomial(m,:),normalizedt);
                end
                X(OCMATCONT.firstordercoordinate,tidx)=coeff1part.'*psival1;
                X(OCMATCONT.zeroordercoordinate,tidx)=coeff0part.'*psival0;
            end
        end
end
