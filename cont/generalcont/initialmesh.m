function coeff=initialmesh(arcmesh,arcval,tmesh,arc)
%INITIALMESH interpolates given values with cubic spline interpolation
%and calculates the corresponding coefficients in the Runge-Kutta basis

global OCMATCONT
domainddata=OCMATCONT.DOMAINDDATA;

diffmesh=diff(tmesh);
% remove zero entries
nummeshintv=numel(diffmesh);
counter=0;

%arcdiffmesh=diff(arcmesh);
arcindex=OCMATCONT.HE.arcindex(arc);
numcols=domainddata(arcindex).numcols;
odecoord=domainddata(arcindex).odecoord;
aecoord=domainddata(arcindex).aecoord;
numeq=domainddata(arcindex).numeq;
numode=domainddata(arcindex).numode;
numae=domainddata(arcindex).numae;
numcolscoord=domainddata(arcindex).numcolscoord;
psival=domainddata(arcindex).psivalequi;
psival0=domainddata(arcindex).psivalequi0;
if ~numae % no algebraic equation (implicit control) and ODEs are of order 1
    equidistmesh=[0:1/(numcols):1].';
    interpolationmesh=tmesh(ones(numcols,1),1:end-1)+equidistmesh(1:numcols,ones(1,nummeshintv)).*diffmesh(ones(numcols,1),:);
    interpolationmesh=[interpolationmesh(:).' tmesh(end)];
    interpolatedval=interp1(arcmesh,arcval(odecoord,:).',interpolationmesh,'spline').';
    A=zeros(numcols+1);
    H=ones(numcols+1);
    A(1:numcols+1,1)=1; % coefficient for y_ijk
    A(:,numcolscoord+1)=reshape(psival(1,numcolscoord,1:numcols+1),numcols,numcols+1).'; % coefficient for z_ijlk
    for jj=0:nummeshintv-1
        H(2:numcols+1,2:numcols+1)=diffmesh(jj+1);
        intervaldynVar=interpolatedval(:,jj*(numcols)+1:(jj+1)*(numcols)+1);
        ((A.*H)\intervaldynVar').';
        coeff([(jj*(numeq)*(numcols+1))+1:(jj+1)*((numeq)*(numcols+1))]+counter,1)=ans(:);
    end
    counter=(jj+1)*((numeq)*(numcols+1))+counter;


else
    % for ODEs and AEs
    equidistmesh=[0:1/(numcols):1].';
    interpolationmesh=tmesh(ones(numcols,1),1:end-1)+equidistmesh(1:numcols,ones(1,nummeshintv)).*diffmesh(ones(numcols,1),:);
    interpolationmesh=[interpolationmesh(:).' tmesh(end)];
    interpolatedval=interp1(arcmesh,arcval(odecoord,:).',interpolationmesh,'spline').';
    equidistmeshae=[0:1/(numcols-1):1].';
    interpolationmeshae=tmesh(ones(numcols-1,1),1:end-1)+equidistmeshae(1:numcols-1,ones(1,nummeshintv)).*diffmesh(ones(numcols-1,1),:);
    interpolationmeshae=[interpolationmeshae(:).' tmesh(end)];
    if numae==1
        interpolatedvalae=interp1(arcmesh,arcval(aecoord,:).',interpolationmeshae,'spline');
    else
        interpolatedvalae=interp1(arcmesh,arcval(aecoord,:).',interpolationmeshae,'spline').';
    end
    Aae=reshape(psival0(1,numcolscoord,1:numcols),numcols,numcols).';
    A=zeros(numcols+1);
    H=ones(numcols+1);
    A(1:numcols+1,1)=1; % coefficient for y_ijk
    A(:,numcolscoord+1)=reshape(psival(1,numcolscoord,1:numcols+1),numcols,numcols+1).'; % coefficient for z_ijlk
    for jj=0:nummeshintv-1
        counter_start=counter+1;
        counter=counter+numode*(numcols+1)+numae*numcols;
        H(2:numcols+1,2:numcols+1)=diffmesh(jj+1);
        intervaldynVar=interpolatedval(odecoord,jj*(numcols)+1:(jj+1)*(numcols)+1);
        Code=((A.*H)\intervaldynVar').';
        intervaldynVarae=interpolatedvalae(1:numae,jj*(numcols-1)+1:(jj+1)*(numcols-1)+1);
        Cae=((Aae)\intervaldynVarae').';
        C=[Code(odecoord,2:numcols+1);Cae];
        coeff(counter_start:counter,1)=[Code(odecoord,1);C(:)];
    end
end


