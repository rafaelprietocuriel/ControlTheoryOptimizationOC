function b=ispdeasymptotic(varargin)

b=cellfun(@isa,varargin,repmat({'pdeasymptotic'},1,nargin),'UniformOutput',1);