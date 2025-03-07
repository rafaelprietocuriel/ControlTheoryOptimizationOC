function b=isocgradlimitset(varargin)

b=cellfun(@isa,varargin,repmat({'ocgradlimitset'},1,nargin),'UniformOutput',1);