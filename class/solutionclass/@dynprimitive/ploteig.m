function varargout=ploteig(dynPrim,varargin)
%
% PLOTEIG plots the eigenvectors

fac=[];
if ~isequilibrium(dynPrim)
    return
end

if nargin>1
    fac=varargin{1};
end
if isempty(fac)
    fac=1;
end

[eigvec,eigval]=reig(dynPrim);
eigval=diag(eigval);
[dum,maxidx]=max(eigval);
[dum,minidx]=min(eigval);
x=dynPrim.octrajectory.y(1:2);
h=zeros(size(eigvec,2),1);
for ii=1:size(eigvec,2)
    if eigval(ii)<0
        clr=[1 0 0];
        totfac=fac;
        linestyle='-';
    elseif abs(eigval(ii))<1e-6
        clr=[0 0 0];
        totfac=0.25*fac;
        linestyle='-';
    else
        clr=[0 0 1];
        if ii==maxidx
            linestyle='-';
        else
            linestyle='--';
        end
        totfac=0.5*fac;
    end
    h(ii)=plot([x(1)-totfac*eigvec(1,ii),x(1)+totfac*eigvec(1,ii)],[x(2)-totfac*eigvec(2,ii),x(2)+totfac*eigvec(2,ii)]);
    set(h(ii),'Color',clr,'LineStyle',linestyle);
end

if nargout>=1
    varargout{1}=h;
end
if nargout>=2
    varargout{2}=eigval;
end
if nargout>=3
    varargout{3}=eigvec;
end