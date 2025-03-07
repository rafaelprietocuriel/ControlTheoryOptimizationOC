function varargout=eig(ppdePrim,varargin)
%
%
if nargout==0
    varargout{1}=eig(full(linearization(ppdePrim)));
else
    [varargout{1:nargout}]=eig(full(linearization(ppdePrim)));
end