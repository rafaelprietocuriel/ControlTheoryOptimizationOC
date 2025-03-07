function b=isocgradasymptotic(varargin)

b=cellfun(@isa,varargin,repmat({'ocgradasymptotic'},1,nargin),'UniformOutput',1);