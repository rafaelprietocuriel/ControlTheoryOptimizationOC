function ocAsy = manifoldpath(ocObj,dynPrim,timelength,varargin)
%
% STABLEPATH calculates the stable path
%
% STABLEPATH(OCOBJ,DYNPRIM,T) calculates the stable path of DYNPRIM with
% time length T using the IVP approach.
%
% STABLEPATH(OCOBJ,DYNPRIM,T,OPT) the options structure OPT provides the
% tolerances for finding eigenvalues of the center manifold (ZeroEigValTol)
% and the distance from the equilibrium (EigTol).  
%
% STABLEPATH(OCOBJ,DYNPRIM,T,OPT,W) if W is set to 'strong' than  the
% strong stable path is calculated (eigenvector corresponding to the
% eigenvalue with maximum absolut real value). Otherwise W is a vector
% denoting the weights for the eigenvectors.

opt=[];
ocAsy=ocasymptotic;
manifoldtype=''; % possible values 's', 'u', 'c', 'sts', 'stu'
weighteigvec=[];
timedirection=1;

if nargin>=4
    opt=varargin{1};
end

if nargin>=5
    manifoldtype=varargin{2};
end

if nargin>=6
    weighteigvec=varargin{3};
end

if isempty(manifoldtype)
    manifoldtype='s';
end

if isempty(opt)
    opt=defaultocoptions;
end

if isempty(ocObj) || isempty(dynPrim)
    return
end

timelength=abs(timelength);
if isempty(timelength) || ~timelength
    return
end

direction=1+(-2)*getocoptions(opt,'OCCONTARG','Backward');
if ~isa(dynPrim,'dynprimitive')
    error('Second input argument must be a dynprimitive');
end

if length(timelength)>1
    error('Incorrect length of input argument timelength');
end
evtol=getocoptions(opt,'OCCONTARG','InitStepWidth');

dynPrim0=dynPrim;

[numseigval,numueigval,numceigval,seigenvector,seigenvalue,ueigenvector,ueigenvalue,ceigenvector,ceigenvalue]=characteristics(dynPrim,opt);

switch lower(manifoldtype(end))
    case 's'
        timedirection=-1;
        num=numseigval{:};
        eigvec=seigenvector{:};
        eigval=seigenvalue{:};
    case 'u'
        num=numueigval{:};
        eigvec=ueigenvector{:};
        eigval=ueigenvalue{:};
    case 'c'
        timedirection=-1;
        num=numceigval{:};
        eigvec=ceigenvector{:};
        eigval=ceigenvalue{:};
    otherwise
end

timelength=timelength*timedirection;

if isempty(eigvec)
    return
end

numeigvec=numel(eigval);
if numel(manifoldtype)>=2
    switch lower(manifoldtype(end-1:end))
        case 'ts'
            weighteigvec=zeros(1,numeigvec);
            [dum idx]=max(abs(real(eigval)));
            %eigvec=eigvec(:,idx);
            %eigval=eigval(idx);
            weighteigvec(idx)=1;
        case 'ws'
            weighteigvec=zeros(1,numeigvec);
            [dum idx]=min(abs(real(eigval)));
            weighteigvec(idx)=1;
    end
elseif isempty(weighteigvec)
    weighteigvec=ones(1,numeigvec);
end

if numel(weighteigvec)~=numeigvec
    error('The vector weighting the eigenvectors has wrong size.')
end
weighteigvec=weighteigvec(:).';
weighteigvec=weighteigvec/norm(weighteigvec);
weighteigvec=repmat(weighteigvec,size(eigvec,1),1)*evtol*direction;


ocTrj=odesolve(ocObj,[0 timelength],dynPrim.y+sum(weighteigvec.*eigvec,2),opt,dynPrim.arcarg,varargin{4:end});
ocAsy=ocasymptotic(ocTrj,dynPrim0);

