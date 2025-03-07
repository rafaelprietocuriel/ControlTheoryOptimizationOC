function varargout=feval(odeObj,varargin)
%
% overloaded function to evaluate function which start with the modelname
% of odeObj.

varargout=cell(1,nargout);
if isempty(odeObj)
    return
end
if nargin==1
    ocmaterror('Not enough input arguments.')
end
if nargout==0
    feval([modelname(odeObj) varargin{1}],varargin{2:end});
    try
        varargout{1}=ans;
    end
else
    [varargout{1:nargout}]=feval([modelname(odeObj) varargin{1}],varargin{2:end});
end