function varargout=loglogcont(varargin)

varargout{1:nargout}=plotcont(varargin{:});
set(gca,'XScale','log','YScale','log')