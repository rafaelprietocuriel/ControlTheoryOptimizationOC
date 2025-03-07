function varargout=eig(pdePrim,varargin)
%
%
if nargout==0
    varargout{1}=eig(full(linearization(pdePrim)));
else
    [varargout{1:nargout}]=eig(full(linearization(pdePrim)));
end