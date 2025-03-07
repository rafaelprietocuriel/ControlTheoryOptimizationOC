function b=isppdeasymptotic(varargin)

b=cellfun(@isa,varargin,repmat({'ppdeasymptotic'},1,nargin),'UniformOutput',1);