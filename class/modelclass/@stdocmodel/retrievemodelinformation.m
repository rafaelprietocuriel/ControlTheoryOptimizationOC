function data=retrievemodelinformation(ocObj,varargin)
%
% RETRIEVEMODELINFORMATION returns a structure with information about the
% ocmat model
% this function is the interface to the structure OCSTRUCT as derived
% from the models initialization file. Its purpose is to allow the same
% commands even if the structure OCSTRUCT is changed. In that case only
% 'retrievemodelinformation' has to be adapted.

if isempty(ocObj)
    data=[];
    return
end
load(ocObj,'modeldata');
data=retrievemodelinformation(ocStruct,varargin{:});