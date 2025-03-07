function b=isdoctrajectory(varargin)

b=cellfun(@isa,varargin,repmat({'doctrajectory'},1,nargin),'UniformOutput',1);
