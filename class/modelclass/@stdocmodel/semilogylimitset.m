function varargout=semilogylimitset(varargin)

varargout{1:nargout}=plotlimitset(varargin{:});
set(gca,'YScale','log')
