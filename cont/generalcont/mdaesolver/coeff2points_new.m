function X=coeff2points_new(tmesh,coeff,flag,varargin)
global OCMATCONT

arcposition=getarcposition(tmesh);

htot=diff(tmesh);
N=diff(arcposition);

%nf=OCMATCONT.firstordercomponentnumber;
%
%N=OCMATCONT.MESHDATA.meshNumber;m=OCMATCONT.CollocationNumber;n=OCMATCONT.componentnumber;nz=OCMATCONT.zeroordercomponentnumber;
%fc=OCMATCONT.firstordercoordinate;zc=OCMATCONT.zeroordercoordinate
switch flag

    case 'general'
        if nargin==3
            return
        end
        t=varargin{1};
        X=zeros(OCMATCONT.componentnumber,length(t));

        ctrcoeff=0;
        for arc=1:OCMATCONT.arcnumber
            [coeff1,coeff0]=splitcoefficient(coeff(ctrcoeff+(1:N(arc)*(OCMATCONT.firstordercomponentnumber+OCMATCONT.componentnumber*OCMATCONT.CollocationNumber))),N(arc)+1,arc);
            ctrcoeff=ctrcoeff+N(arc)*(OCMATCONT.firstordercomponentnumber+OCMATCONT.componentnumber*OCMATCONT.CollocationNumber);

            arctmesh=tmesh(arcposition(1,arc):arcposition(2,arc));
            h=htot(arcposition(1,arc):arcposition(2,arc)-1);
            for ii=1:OCMATCONT.MESHDATA.meshNumber(arc)-1
                if ii<OCMATCONT.MESHDATA.meshNumber(arc)-1
                    % exclude the right side of the interval
                    tidx=find(arctmesh(ii)<=t & t<arctmesh(ii+1));
                else
                    % include the right side of the interval
                    tidx=find(arctmesh(ii)<=t & t<=arctmesh(ii+1));
                end
                tidx=tidx+arcposition(1,arc)-1;
                normalizedt=(t(tidx)-tmesh(ii))/h(ii);
                coeff1part=coeff1(:,(ii-1)*OCMATCONT.firstordercomponentnumber(arc)+(1:OCMATCONT.firstordercomponentnumber(arc)));
                coeff0part=coeff0(:,(ii-1)*OCMATCONT.zeroordercomponentnumber(arc)+(1:OCMATCONT.zeroordercomponentnumber(arc)));
                psival1=ones(OCMATCONT.CollocationNumber+1,length(tidx));
                psival0=zeros(OCMATCONT.CollocationNumber,length(tidx));
                for jj=1:OCMATCONT.CollocationNumber
                    psival1(jj+1,:)=h(ii)*polyval(OCMATCONT.FirstOrderCollocationPolynomial(jj,:),normalizedt);
                    psival0(jj,:)=polyval(OCMATCONT.ZeroOrderCollocationPolynomial(jj,:),normalizedt);
                end
                X(OCMATCONT.firstordercoordinate,tidx)=coeff1part.'*psival1;
                X(OCMATCONT.zeroordercoordinate,tidx)=coeff0part.'*psival0;
            end
        end
        
    case 'grid'
        X=zeros(OCMATCONT.componentnumber,sum(OCMATCONT.MESHDATA.meshNumber));

        ctr=0;
        ctrcoeff=0;
        for arc=1:OCMATCONT.arcnumber
            [coeff1,coeff0]=splitcoefficient(coeff(ctrcoeff+(1:N(arc)*(OCMATCONT.firstordercomponentnumber+OCMATCONT.componentnumber*OCMATCONT.CollocationNumber))),N(arc)+1,arc);
            
            H=myrepmat(htot(arcposition(1,arc):arcposition(2,arc)-1),OCMATCONT.firstordercomponentnumber,1,1,N(arc));
            H=myrepmat(H(:).',OCMATCONT.CollocationNumber,1,1,N(arc)*OCMATCONT.firstordercomponentnumber);
            coeff1(2:OCMATCONT.CollocationNumber+1,:)=coeff1(2:OCMATCONT.CollocationNumber+1,:).*H;
            X1=reshape(OCMATCONT.FirstOrderCollocationValues(:,1).'*coeff1,OCMATCONT.firstordercomponentnumber,[]);
            X0=reshape(OCMATCONT.ZeroOrderCollocationValues(:,1).'*coeff0,OCMATCONT.zeroordercomponentnumber,[]);
            X1T=OCMATCONT.FirstOrderCollocationValues(:,OCMATCONT.CollocationNumber+2).'*coeff1(:,end-OCMATCONT.firstordercomponentnumber+1:end);
            X0T=OCMATCONT.ZeroOrderCollocationValues(:,OCMATCONT.CollocationNumber+2).'*coeff0(:,end-OCMATCONT.zeroordercomponentnumber+1:end);
            X(OCMATCONT.firstordercoordinate,ctr+(1:OCMATCONT.MESHDATA.meshNumber(arc)))=[X1 X1T(:)];
            X(OCMATCONT.zeroordercoordinate,ctr+(1:OCMATCONT.MESHDATA.meshNumber(arc)))=[X0 X0T(:)];
            if arc<OCMATCONT.arcnumber
                ctrcoeff=ctrcoeff+N(arc)*(OCMATCONT.firstordercomponentnumber+OCMATCONT.componentnumber*OCMATCONT.CollocationNumber);
                ctr=ctr+((OCMATCONT.MESHDATA.meshNumber(arc)-1)*(OCMATCONT.CollocationNumber+1)+1);
            end
        end

        
    case 'collocationgrid'
        % the components are computed als the polynomial values evaluated
        % at 0 rho1 ... rho_m
        X=zeros(OCMATCONT.componentnumber,sum(N)*(OCMATCONT.CollocationNumber+1)+OCMATCONT.arcnumber);

        ctrcoeff=0;
        ctr=0;
        for arc=1:OCMATCONT.arcnumber
            Narc=N(arc)+1;
            [coeff1,coeff0]=splitcoefficient(coeff(ctrcoeff+(1:N(arc)*(OCMATCONT.firstordercomponentnumber+OCMATCONT.componentnumber*OCMATCONT.CollocationNumber))),Narc,arc);
            
            H=myrepmat(htot(arcposition(1,arc):arcposition(2,arc)-1),OCMATCONT.firstordercomponentnumber,1,1,N(arc));
            H=myrepmat(H(:).',OCMATCONT.CollocationNumber,1,1,N(arc)*OCMATCONT.firstordercomponentnumber);
            coeff1(2:OCMATCONT.CollocationNumber+1,:)=coeff1(2:OCMATCONT.CollocationNumber+1,:).*H;
            X1=reshape(permute(reshape(OCMATCONT.FirstOrderCollocationValues(:,1:OCMATCONT.CollocationNumber+1).'*coeff1,OCMATCONT.CollocationNumber+1,OCMATCONT.firstordercomponentnumber,[]),[2 1 3]),OCMATCONT.firstordercomponentnumber,[]);
            X0=reshape(permute(reshape(OCMATCONT.ZeroOrderCollocationValues(:,1:OCMATCONT.CollocationNumber+1).'*coeff0,OCMATCONT.CollocationNumber+1,OCMATCONT.zeroordercomponentnumber,[]),[2 1 3]),OCMATCONT.zeroordercomponentnumber,[]);
            X1T=OCMATCONT.FirstOrderCollocationValues(:,OCMATCONT.CollocationNumber+2).'*coeff1(:,end-OCMATCONT.firstordercomponentnumber+1:end);
            X0T=OCMATCONT.ZeroOrderCollocationValues(:,OCMATCONT.CollocationNumber+2).'*coeff0(:,end-OCMATCONT.zeroordercomponentnumber+1:end);
            X(OCMATCONT.firstordercoordinate,ctr+(1:((OCMATCONT.MESHDATA.meshNumber(arc)-1)*(OCMATCONT.CollocationNumber+1)+1)))=[X1 X1T(:)];
            X(OCMATCONT.zeroordercoordinate,ctr+(1:((OCMATCONT.MESHDATA.meshNumber(arc)-1)*(OCMATCONT.CollocationNumber+1)+1)))=[X0 X0T(:)];
            if arc<OCMATCONT.arcnumber
                ctrcoeff=ctrcoeff+N(arc)*(OCMATCONT.firstordercomponentnumber+OCMATCONT.componentnumber*OCMATCONT.CollocationNumber);
                ctr=ctr+((OCMATCONT.MESHDATA.meshNumber(arc)-1)*(OCMATCONT.CollocationNumber+1)+1);
            end
        end

end