function varargout=feval(dgObj,varargin)
%
% overloaded function to evaluate function which start with the modelname
% of dgObj.

varargout=cell(1,nargout);
if isempty(dgObj)
    return
end
if nargin==1
    ocmaterror('Not enough input arguments.')
end
if nargout==0
    feval([modelname(dgObj) varargin{1}],varargin{2:end});
    try
        varargout{1}=ans;
    end
else
    [varargout{1:nargout}]=feval([modelname(dgObj) varargin{1}],varargin{2:end});
end