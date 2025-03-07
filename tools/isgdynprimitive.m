function b=isgdynprimitive(varargin)

%b=cellfun(@isa,varargin,repmat({'dynprimitive'},1,nargin),'UniformOutput',1);
b=strcmp(cellfun(@class,varargin,'UniformOutput',0),'gdynprimitive');