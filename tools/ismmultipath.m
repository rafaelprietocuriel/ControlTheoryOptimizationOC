function b=ismmultipath(varargin)

b=cellfun(@isa,varargin,repmat({'mmultipath'},1,nargin),'UniformOutput',1);
