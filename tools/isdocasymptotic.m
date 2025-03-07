function b=isdocasymptotic(varargin)

b=cellfun(@isa,varargin,repmat({'docasymptotic'},1,nargin),'UniformOutput',1);