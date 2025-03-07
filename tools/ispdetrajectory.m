function b=ispdetrajectory(varargin)

b=cellfun(@isa,varargin,repmat({'pdetrajectory'},1,nargin),'UniformOutput',1);
