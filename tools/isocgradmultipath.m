function b=isocgradmultipath(varargin)

b=cellfun(@isa,varargin,repmat({'ocgradmultipath'},1,nargin),'UniformOutput',1);
