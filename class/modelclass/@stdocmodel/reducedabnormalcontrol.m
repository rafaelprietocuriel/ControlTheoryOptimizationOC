function varargout=reducedabnormalcontrol(ocObj,X,arcid)
%
if isempty(ocObj)
    varargout{1}=[];
    return
end
symbolic=false;
if nargin==1 || isempty(X)
    symbolic=true;
end
if nargin<3
    arcid=[];
end
if nargin<4
    repflag=[];
end
if isempty(arcid)
    arcid=0;
end
if symbolic
    varargout{1}=feval(ocObj,'ReducedAbnormalSymbolicOptimalControl',arcid);
else
    varargout{1}=feval(ocObj,'ReducedAbnormalOptimalControl',arcid,repflag);
end