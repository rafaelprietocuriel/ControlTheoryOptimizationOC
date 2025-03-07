function varargout=reig(dynPrim,flag,varargin)
%
% EIG find eigenvalues and eigenvectors.
%
% D = EIG(DYNPRIM) returns the eigenvalues for the Jacobian/Monodromy
% matrix of the dynprimitve object DYNPRIM.    
%
% [D,V] = EIG(DYNPRIM) produces matrices of eigenvalues (D) and eigenvectors
% (V) (for more information help eig).  

if nargin==1
    flag=1; % returns Jacobian calculated for the state costate coordinates
end
L=jacobian(dynPrim,flag);
if nargout==0
    varargout{1}=reig(L{1},varargin{:});
else
    [varargout{1:nargout}]=reig(L{1},varargin{:});
end