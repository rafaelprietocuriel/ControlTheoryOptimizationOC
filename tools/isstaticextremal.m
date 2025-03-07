function b=isstaticextremal(varargin)

b=cellfun(@isa,varargin,repmat({'staticextremal'},1,nargin),'UniformOutput',1);
