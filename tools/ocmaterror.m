function ocmaterror(varargin)
%
%

if nargin>1
    error(['OCMatError:' varargin{1}],varargin{2:end})
else
    error('OCMatError:Standard',varargin{:})
end