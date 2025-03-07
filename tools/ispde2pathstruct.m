function b=ispde2pathstruct(varargin)

b=cellfun(@isfield,varargin,repmat({'lss'},1,nargin),'UniformOutput',1);
