function b=ismapprimitive(varargin)

b=cellfun(@isa,varargin,repmat({'mapprimitive'},1,nargin),'UniformOutput',1);
