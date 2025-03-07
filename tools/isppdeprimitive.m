function b=isppdeprimitive(varargin)

b=cellfun(@isa,varargin,repmat({'ppdeprimitive'},1,nargin),'UniformOutput',1) && ~cellfun(@isa,varargin,repmat({'ppdeasymptotic'},1,nargin),'UniformOutput',1);
