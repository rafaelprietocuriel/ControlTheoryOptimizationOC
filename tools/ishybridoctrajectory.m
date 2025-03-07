function b=ishybridoctrajectory(varargin)

b=cellfun(@isa,varargin,repmat({'hybridoctrajectory'},1,nargin),'UniformOutput',1);