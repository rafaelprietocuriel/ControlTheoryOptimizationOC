function b=isocasymptotic(varargin)

b=cellfun(@isa,varargin,repmat({'ocasymptotic'},1,nargin),'UniformOutput',1);