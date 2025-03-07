function varargout=semilogxlimitset(varargin)

varargout{1:nargout}=plotlimitset(varargin{:});
set(gca,'XScale','log')
