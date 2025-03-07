function b=ispdeprimitive(varargin)

b=cellfun(@isa,varargin,repmat({'pdeprimitive'},1,nargin),'UniformOutput',1) && ~cellfun(@isa,varargin,repmat({'pdeasymptotic'},1,nargin),'UniformOutput',1);
