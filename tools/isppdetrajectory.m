function b=isppdetrajectory(varargin)

b=cellfun(@isa,varargin,repmat({'ppdetrajectory'},1,nargin),'UniformOutput',1);
