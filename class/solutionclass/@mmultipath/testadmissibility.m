function [b infoS newarcpos violarcarg]=testadmissibility(ocMP,ocObj,varargin)

opt=[];
conttype='';
if nargin>=3
    opt=varargin{1};
end
if nargin>=4
    conttype=varargin{2};
end
if isempty(opt)
    opt=defaultocoptions;
end
if isempty(ocMP)
    b=[];
    infoS=[];
    newarcpos=[];
    violarcarg=[];
    return
end

for ii=1:multiplicity(ocMP)
    [b(ii) infoS(ii) newarcpos{ii} violarcarg{ii}]=testadmissibility(ocMP.solutionclass{ii},ocObj,opt,conttype);
end