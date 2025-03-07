function b=ismodelcomposite(varargin)

b=cellfun(@isa,varargin,repmat({'modelcomposite'},1,nargin),'UniformOutput',1);
