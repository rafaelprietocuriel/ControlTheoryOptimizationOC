function varargout=reig(dynPrim,varargin)
%
% EIG find eigenvalues and eigenvectors.
%
% D = EIG(DYNPRIM) returns the eigenvalues for the Jacobian/Monodromy
% matrix of the dynprimitve object DYNPRIM.    
%
% [D,V] = EIG(DYNPRIM) produces matrices of eigenvalues (D) and eigenvectors
% (V) (for more information help eig).  

L=linearization(dynPrim);
if nargout==0
    varargout{1}=reig(L,varargin{:});
else
    [varargout{1:nargout}]=reig(L,varargin{:});
end