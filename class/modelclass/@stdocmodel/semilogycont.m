function varargout=semilogycont(varargin)

varargout{1:nargout}=plotcont(varargin{:});
set(gca,'YScale','log')

