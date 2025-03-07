function X=coeff2derivatives(tcolmesh,coeff,flag,varargin)
global OCMATCONT

switch flag
        
    case 'general'
        % coefficients based on data stored in MESHDATA
        if nargin==2
            return
        end
        t=varargin{1};
        tmesh=tcolmesh(OCMATCONT.MESHDATA.tmeshidx);
        coeff1=coeff(OCMATCONT.MESHDATA.firstorderidx(2:end,:));
        X=zeros(OCMATCONT.componentnumber,length(t));
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
                psival0=zeros(OCMATCONT.CollocationNumber,length(tidx));
                for m=1:OCMATCONT.CollocationNumber
                    psival0(m,:)=polyval(OCMATCONT.ZeroOrderCollocationPolynomial(m,:),normalizedt);
                end
                X(OCMATCONT.firstordercoordinate,tidx)=coeff1part.'*psival0;
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
            X=zeros(OCMATCONT.componentnumber,length(t));
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
                    psival1(m+1,:)=OCMATCONT.MESHDATA1_2.h(jj)*polyval(OCMATCONT.FirstOrderCollocationPolynomial(m,:),normalizedt);
                    psival0(m,:)=polyval(OCMATCONT.ZeroOrderCollocationPolynomial(m,:),normalizedt);
                end
                X(OCMATCONT.firstordercoordinate,tidx)=coeff1part.'*psival1;
                X(OCMATCONT.zeroordercoordinate,tidx)=coeff0part.'*psival0;
            end
        end
end