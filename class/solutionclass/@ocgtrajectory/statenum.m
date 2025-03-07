function out=statenum(ocgTrj,varargin)
%
%
ocObj=loadmodel(ocgTrj);
if isempty(ocObj)
    out=[];
    return
end
out=statenum(ocObj);