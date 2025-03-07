function varargout=semilogxcont(varargin)

varargout{1:nargout}=plotcont(varargin{:});
set(gca,'XScale','log')

