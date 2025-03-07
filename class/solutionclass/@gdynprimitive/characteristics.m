function [numseigval,numueigval,numceigval,seigenvector,seigenvalue,ueigenvector,ueigenvalue,ceigenvector,ceigenvalue]=characteristics(varargin)
%
% CHARACTERISTICS returns the eigenvalues and corresponding eigenvectors
%
% CHARACTERISTICS(DYNPRIM)
%
% CHARACTERISTICS(DYNPRIM,OPT)
%
% CHARACTERISTICS(DYNPRIM1,...,DYNPRIMN,OPT)
%
% [SEVAL,UEVAL,CEVAL,SEVEC,UEVEC,CEVEC]=CHARACTERISTICS(...)

opt=[];
numseigval{1}=[];
numueigval{1}=[];
numceigval{1}=[]; % number of eigenvalues with zero real value
seigenvalue{1}=[];
ueigenvalue{1}=[];
ceigenvalue{1}=[];
seigenvector{1}=[];
ueigenvector{1}=[];
ceigenvector{1}=[];

if nargin>=2
    if isstruct(varargin{end})
        opt=varargin{end};
        varargin(end)=[];
    end
end

if isempty(opt)
    opt=defaultocoptions;
end

tol=getocoptions(opt,'GENERAL','ZeroDeviationTolerance');

numdynprim=length(varargin);
for ii=1:numdynprim
    try
        [eigvec eigval]=eig(varargin{ii});
    catch
        eigval=[];
    end
    if ~isempty(eigval)
        n=length(eigval);
        eigval=diag(eigval);

        if isequilibrium(varargin{ii})
            ff=find(abs(real(eigval))<tol);
            eigval(ff)=0;
            
            ff=find(real(eigval)<0);
            numseigval{ii}=length(ff);
            seigenvalue{ii}=eigval(ff);
            seigenvector{ii}=eigvec(:,ff);

            ff=find(real(eigval)>0);
            numueigval{ii}=length(ff);
            ueigenvalue{ii}=eigval(ff);
            ueigenvector{ii}=eigvec(:,ff);

            ff=find(real(eigval)==0);
            numceigval{ii}=length(ff);
            ceigenvalue{ii}=eigval(ff);
            ceigenvector{ii}=eigvec(:,ff);

        elseif  islimitcycle(varargin{ii})
            ff=find(abs(real(eigval)-1)<tol);
            eigval(ff)=1;

            ff=find(real(eigval)<1);
            numseigval{ii}=length(ff);
            seigenvalue{ii}=eigval(ff);
            seigenvector{ii}=eigvec(:,ff);

            ff=find(real(eigval)>1);
            numueigval{ii}=length(ff);
            ueigenvalue{ii}=eigval(ff);
            ueigenvector{ii}=eigvec(:,ff);

            ff=find(real(eigval)==1);
            numceigval{ii}=length(ff);
            ceigenvalue{ii}=eigval(ff);
            ceigenvector{ii}=eigvec(:,ff);
        end
    else
        numseigval{ii}=0;
        numueigval{ii}=0;
        numceigval{ii}=0;
        ceigenvalue{ii}=[];
        ceigenvector{ii}=[];
        ueigenvalue{ii}=[];
        ueigenvector{ii}=[];
        seigenvalue{ii}=[];
        seigenvector{ii}=[];
    end
end