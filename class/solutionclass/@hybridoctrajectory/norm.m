function out=norm(hocTrj,varargin)

if nargin==1
    flag='t'; % total
else
    flag=varargin{1};
end

if isfield(hocTrj.solverinfo,'objectivevaluecoord')
    hocTrj.y(hocTrj.solverinfo.objectivevaluecoord,:)=[];
end


arcp=arcposition(hocTrj);
jumparg=jumpargument(hocTrj);
t=time(hocTrj,1);
t=t([1 1:end end]);
xm=hocTrj.y(:,1);
dnorm=0;
cnorm=0;
for ii=1:length(jumparg)
    if ii<length(jumparg)
        xr=hocTrj.y(:,arcp(1,ii));
        dx=sqrt(sum(hocTrj.y(:,arcp(1,ii):arcp(2,ii)).^2,1));
        dt=diff(t(arcp(1,ii):arcp(2,ii)));
        cnorm=cnorm+sum((dx(1:end-1)+dx(2:end))/2.*dt);
        dnorm=dnorm+sqrt(sum((xr-xm).^2));
        xm=hocTrj.y(:,arcp(2,ii));
    else
        xr=hocTrj.y(:,end);
        dnorm=dnorm+sqrt(sum((xr-xm).^2));
    end
end

switch flag
    case 't'
        out=dnorm+cnorm;
    case 'c'
        out=cnorm;
    case 'd'
        out=dnorm;
    otherwise
        out=[];
end