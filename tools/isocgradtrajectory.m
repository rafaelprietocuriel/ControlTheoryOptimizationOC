function b=isocgradtrajectory(varargin)

b=cellfun(@isa,varargin,repmat({'ocgradtrajectory'},1,nargin),'UniformOutput',1);
