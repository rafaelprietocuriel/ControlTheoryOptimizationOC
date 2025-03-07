function varargout=reducedabnormalcanonicalsystem(ocObj,arcid,repflag)
%
if isempty(ocObj)
    varargout{1}=[];
    return
end
if nargin<2
    arcid=[];
end
if nargin<3
    repflag=[];
end
if isempty(arcid)
    arcid=0;
end
if isempty(repflag)
    repflag=0;
end
varargout{1}=feval(ocObj,'ReducedAbnormalSymbolicCanonicalSystem',arcid,repflag);
