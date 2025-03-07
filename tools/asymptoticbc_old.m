function [F,numstable,numunstable,numcenter,infoStruct]=asymptoticbc(J,varargin)
%
% ASYMPTOTICBC returns the matrix for the asymptotic boundary condition
%
% see W.-J. Beyn, T. Pampel, W.Semmler, 2001. "Dynamic optimization and
% Skiba sets in economic examples," Computing in Economics and Finance 2001
% 29, Society for Computational Economics. p. 4 (calculation of V)

F=[];
numstable=[];
numunstable=[];
numcenter=[];
infoStruct=[];
limitsettype='';
stabilitytype='';
centertol=[];
if isempty(J)
    return
end
if nargin>=2
    stabilitytype=varargin{1};
end
if nargin>=3
    limitsettype=varargin{2};
end
if nargin>=4
    centertol=varargin{3};
end
if isempty(stabilitytype)
    stabilitytype='s';
end
if isempty(limitsettype)
    limitsettype='c'; % (c)ontinuous [default] or (d)iscrete
end
if isempty(centertol)
    centertol=1e-8;
end

[U,D,eigval]=realeig(J);
reigval=real(eigval);
% uses QR factorization for the calculation of the the matrix Q, that is
% orthogonal to the according eigenspace.
%  idx ... the position of the according eigenvectors in U
%   [Q,dum]=qr(U(:,idx));
%Alternatively the SVD decomposition can be used.
%   [Q,dum dum]=svd(U(:,idx));
% B=Q(:,1:length(idx)) ... orthogonal basis of the acccording eigenspace
% F=Q(:,length(idx)+1:end) ... matrix orthogonal to the according
% eigenspace
%
%Another alternative is the null function yielding the nullspace of a
%matrix A, null(A)

n=size(J,1);
switch (stabilitytype)
    case 's'
        if strcmp(limitsettype,'c')
            idx=find(reigval<-centertol);
            [Q,dum]=qr(U(:,idx));
            %F=null(U(:,idx).');
            F=Q(:,length(idx)+1:end);
            numcenter=numel(find(abs(reigval)<centertol));
            numstable=numel(idx);
            numunstable=n-numstable-numcenter;
        elseif strcmp(limitsettype,'d')
            idx=find(abs(reigval)<1-centertol);
            %idx=(abs(reigval)<1-centertol);
            %[Q,TS]=ordschur(Qt,RUt,idx);

            [Q,dum]=qr(U(:,idx));
            %F=null(U(:,idx).');
            F=Q(:,length(idx)+1:end);
            numstable=numel(idx);
            numcenter=numel(find(abs(reigval-1)<centertol));
            numunstable=n-numstable-numcenter;
        end
    case 'S'
        if strcmp(limitsettype,'c')
            idx=find(reigval<0);
            [Q,dum]=qr(U(:,idx));
            F=Q(:,length(idx)+1:end);
            numcenter=numel(find(reigval<centertol & reigval>=0));
            numstable=numel(idx);
            numunstable=n-numstable-numcenter;
        elseif strcmp(limitsettype,'d')
            idx=find(abs(reigval)<1-centertol);
            [Q,dum]=qr(U(:,idx));
            F=Q(:,length(idx)+1:end);
            numstable=numel(idx);
            numcenter=numel(find(abs(reigval-1)<centertol));
            numunstable=n-numstable-numcenter;
        end
    case 'u'
        if strcmp(limitsettype,'c')
            idx=find(reigval>centertol);
            [Q,dum]=qr(U(:,idx));
            F=Q(:,length(idx)+1:end);
            numunstable=numel(idx);
            numcenter=numel(find(abs(reigval)<centertol));
            numstable=n-numunstable-numcenter;
        elseif strcmp(limitsettype,'d')
            idx=find(abs(reigval)>1+centertol);
            [Q,dum]=qr(U(:,idx));
            F=Q(:,length(idx)+1:end);
            numstable=numel(idx);
            numcenter=numel(find(abs(reigval-1)<centertol));
            numunstable=n-numstable-numcenter;
        end
    case 's'
        if strcmp(limitsettype,'c')
            idx=find(reigval<-centertol);
            [Q,dum]=qr(U(:,idx));
            %F=null(U(:,idx).');
            F=Q(:,length(idx)+1:end);
            numcenter=numel(find(abs(reigval)<centertol));
            numstable=numel(idx);
            numunstable=n-numstable-numcenter;
        elseif strcmp(limitsettype,'d')
            idx=find(abs(reigval)<1-centertol);
            [Q,dum]=qr(U(:,idx));
            %F=null(U(:,idx).');
            F=Q(:,length(idx)+1:end);
            numstable=numel(idx);
            numcenter=numel(find(abs(reigval-1)<centertol));
            numunstable=n-numstable-numcenter;
        end
    case 'os'
        if strcmp(limitsettype,'c')
            Q=[];
            idx=find(reigval<-centertol);
            F=U(:,idx);
            numcenter=numel(idx);
            numunstable=numel(find(reigval>centertol));
            numstable=n-numunstable-numcenter;
        elseif strcmp(limitsettype,'d')
            idx=find(abs(abs(reigval)-1)<centertol);
            [Q,dum]=qr(U(:,idx));
            F=Q(:,length(idx)+1:end);
            numcenter=numel(idx);
            numunstable=numel(find(abs(reigval)>1+centertol));
            numstable=n-numunstable-numcenter;
        else
        end
    case {'sc','cs'}
        if strcmp(limitsettype,'c')
            idx=find(reigval<centertol);
            [Q,dum]=qr(U(:,idx));
            F=Q(:,length(idx)+1:end);
            numcenter=numel(find(abs(reigval)<centertol));
            numstable=numel(idx)-numcenter;
            numunstable=n-numstable-numcenter;
        elseif strcmp(limitsettype,'d')
            idx=find(abs(abs(reigval))<1+centertol);
            [Q,dum]=qr(U(:,idx));
            F=Q(:,length(idx)+1:end);
            numcenter=numel(find(abs(abs(reigval)-1)<centertol));
            numunstable=numel(idx)-numcenter;;
            numstable=n-numunstable-numcenter;
        end
    case {'uc','cu'}
        if strcmp(limitsettype,'c')
            idx=find(reigval>-centertol);
            [Q,dum]=qr(U(:,idx));
            F=Q(:,length(idx)+1:end);
            numcenter=numel(find(abs(reigval)<centertol));
            numunstable=numel(idx)-numcenter;
            numstable=n-numunstable-numcenter;
        elseif strcmp(limitsettype,'d')
            eigval=diag(T);
            idx=find(abs(abs(eigval)-1)>centertol & abs(eigval)>1);
            F=U(:,idx);
        end
    case {'sts'}
        if strcmp(limitsettype,'c')
            [minval idx]=min(reigval);
            if minval>=0
                idx=[];
            end
            [Q,dum]=qr(U(:,idx));
            F=Q(:,length(idx)+1:end);
            numcenter=numel(find(abs(reigval)<centertol));
            numstable=numel(idx)-numcenter;
            numunstable=n-numstable-numcenter;
        end
    case {'stu'}
        if strcmp(limitsettype,'c')
            [minval idx]=max(reigval);
            if minval<0
                idx=[];
            end
            [Q,dum]=qr(U(:,idx));
            F=Q(:,length(idx)+1:end);
            numcenter=numel(find(abs(reigval)<centertol));
            numstable=numel(idx)-numcenter;
            numunstable=n-numstable-numcenter;
        end
    case {'ws'}
        if strcmp(limitsettype,'c')
            reigval(reigval>=0)=-inf;
            [maxval idx]=max(reigval);
            if isinf(maxval)
                idx=[];
            end
            [Q,dum]=qr(U(:,idx));
            F=Q(:,length(idx)+1:end);
            numcenter=numel(find(abs(reigval)<centertol));
            numstable=numel(idx)-numcenter;
            numunstable=n-numstable-numcenter;
        end
    case {'wu'}
        if strcmp(limitsettype,'c')
            reigval(reigval<=0)=inf;
            [maxval,idx]=min(reigval);
            if isinf(maxval)
                idx=[];
            end
            [Q,dum]=qr(U(:,idx));
            F=Q(:,length(idx)+1:end);
            numcenter=numel(find(abs(reigval)<centertol));
            numstable=numel(idx)-numcenter;
            numunstable=n-numstable-numcenter;
        end
end
infoStruct.eigval=eigval;
infoStruct.eigvec=U;
infoStruct.index=idx;
infoStruct.Q=Q; % orthogonal basis of the according eigenspace
