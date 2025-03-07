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
method='';
Q=[];
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
if nargin>=5
    method=varargin{4};
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
if isempty(method)
    method='schur';
end

switch method
    case 'schur'
        [Qt,RUt]=schur(J);
        reigval=diag(RUt);
    otherwise
        [U,D,reigval]=realeig(J);
end
%uses QR factorization for the calculation of the the matrix Q, that is
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
            idx=reigval<-centertol;

            numstable=numel(find(idx));
            numcenter=numel(find(abs(reigval)<centertol));
            numunstable=n-numstable-numcenter;
        elseif strcmp(limitsettype,'d')
            idx=(abs(reigval)<1+centertol);

            numstable=numel(find(idx));
            numcenter=numel(find(abs(reigval-1)<centertol));
            numunstable=n-numstable-numcenter;
        end
    case 'u'
        if strcmp(limitsettype,'c')
            idx=reigval>centertol;

            numunstable=sum(idx);
            numcenter=numel(find(abs(reigval)<centertol));
            numstable=n-numunstable-numcenter;
        elseif strcmp(limitsettype,'d')
            idx=abs(reigval)>1+centertol;

            numunstable=sum(idx);
            numcenter=numel(find(abs(reigval-1)<centertol));
            numstable=n-numunstable-numcenter;
        end
    case {'sc','cs'}
        if strcmp(limitsettype,'c')
            idx=reigval<centertol;

            numcenter=numel(find(abs(reigval)<centertol));
            numstable=numel(find(idx))-numcenter;
            numunstable=n-numstable-numcenter;
        elseif strcmp(limitsettype,'d')
            idx=abs(reigval)<1+centertol;
            numstable=numel(find((abs(reigval)<1-centertol)));
            numcenter=numel(find(abs(reigval-1)<centertol));
            numunstable=n-numstable-numcenter;
        end
    case 'c'
        if strcmp(limitsettype,'c')
            idx=abs(reigval)<centertol;

            numcenter=numel(find(abs(reigval)<centertol));
            numstable=numel(find(idx))-numcenter;
            numunstable=n-numstable-numcenter;
        elseif strcmp(limitsettype,'d')
            idx=abs(reigval)<1+centertol;
            numstable=numel(find((abs(reigval)<1-centertol)));
            numcenter=numel(find(abs(reigval-1)<centertol));
            numunstable=n-numstable-numcenter;
        end
    case {'ws'}
        if strcmp(limitsettype,'c')
            reigval(reigval>=0)=-inf;
            maxval=max(reigval);
            idx=reigval==maxval;
            if isinf(maxval)
                idx=[];
            end
            numcenter=numel(find(abs(reigval)<centertol));
            numstable=numel(find(idx))-numcenter;
            numunstable=n-numstable-numcenter;
        end
    case {'sts'}
        if strcmp(limitsettype,'c')
            reigval(reigval>=0)=+inf;
            minval=min(reigval);
            idx=reigval==minval;
            if isinf(minval)
                idx=[];
            end
            numcenter=numel(find(abs(reigval)<centertol));
            numstable=numel(find(idx))-numcenter;
            numunstable=n-numstable-numcenter;
        end
    case {'stu'}
        if strcmp(limitsettype,'c')
            maxval=max(reigval);
            idx=reigval==maxval;
            numcenter=numel(find(abs(reigval)<centertol));
            numstable=numel(find(idx))-numcenter;
            numunstable=n-numstable-numcenter;
        end
end
switch method
    case 'schur'
        [Qt,RUt]=schur(J);
        [Q,TS]=ordschur(Qt,RUt,idx);
        F=Q(:,sum(idx)+1:end);
    case 'qr'
        [Q,dum]=qr(U(:,idx));
        F=Q(:,sum(idx)+1:end);
end

infoStruct.eigval=reigval;
infoStruct.index=idx;
infoStruct.Q=Q; % orthogonal basis of the according eigenspace
