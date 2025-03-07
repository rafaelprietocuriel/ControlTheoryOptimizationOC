function [b,infoS,newarcpos,violarcarg]=testadmissibility(ocAsym,ocObj,varargin)

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
if isempty(conttype)
    conttype='extremal2ep';
end
if isempty(ocAsym)
    b=[];
    infoS=[];
    newarcpos=[];
    violarcarg=[];
    return
end

[b infoS newarcpos violarcarg]=testadmissibility(octrajectory(ocAsym),ocObj,opt,conttype);