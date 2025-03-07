function b=isdocmultipath(varargin)

b=cellfun(@isa,varargin,repmat({'docmultipath'},1,nargin),'UniformOutput',1);
