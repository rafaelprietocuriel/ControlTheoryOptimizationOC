function b=ishybridocmultipath(varargin)

b=cellfun(@isa,varargin,repmat({'hybridocmultipath'},1,nargin),'UniformOutput',1);
