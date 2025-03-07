function [X,par,ocEP,cstr,lm]=testequilibrium(ocObj,matRes,idx,opt)
if nargin==3
    opt=[];
end
if isempty(opt)
    opt=defaultocoptions;
end
matRes.ContinuationSolution.modelparameter=parametervalue(ocObj);
[X,par,idx,arcid]=matcontsolution2ep(matRes,[],idx);

ocObj=changeparametervalue(ocObj,par);

ocEP=calcep(ocObj,X,arcid,opt);
if isempty(ocEP)
    cstr=[];
    lm=[];
    return
end
ocEP=ocEP{1};

if strcmp(getocoptions(opt,'INIT','MessageDisplay'),'on')
    [cstr,lm]=showdata(ocEP);
else
    if isa(ocEP,'gdynprimitive')
        lm=lagrangemultiplier(ocEP);
        lm=lm{1};
        cstr=constraint(ocEP);
        cstr=cstr{1};
    else
        lm=lagrangemultiplier(ocObj,ocEP);
        cstr=constraint(ocObj,ocEP);
    end
end