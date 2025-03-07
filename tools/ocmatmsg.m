function ocmatmsg(varargin)
%
%

if nargin == 0
  error('ocmatmsg needs at least one parameter.');
end

if nargin == 1
  fprintf(varargin{1});
else
  fprintf(varargin{1}, varargin{2:nargin});
end
