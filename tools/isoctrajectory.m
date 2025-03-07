function b=isoctrajectory(varargin)

b=cellfun(@isa,varargin,repmat({'octrajectory'},1,nargin),'UniformOutput',1);
