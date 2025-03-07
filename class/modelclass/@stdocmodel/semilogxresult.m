function varargout=semilogxresult(varargin)

varargout{1:nargout}=plotresult(varargin{:});
set(gca,'XScale','log')

