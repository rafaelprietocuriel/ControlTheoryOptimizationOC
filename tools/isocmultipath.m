function b=isocmultipath(varargin)

b=cellfun(@isa,varargin,repmat({'ocmultipath'},1,nargin),'UniformOutput',1);
