function b=isoccomposite(varargin)

b=cellfun(@isa,varargin,repmat({'occomposite'},1,nargin),'UniformOutput',1);
