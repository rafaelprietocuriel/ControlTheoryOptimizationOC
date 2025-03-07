function fidwrite=remark(odeObj,varargin)
%
% REMARK

%
% REMARK creates or opens a text file for remarks.
%
% REMARK(ODEOBJ) creates a file with the basic filename of the model ODEOBJ,
% see ocmat/tools/remark for further information.

folder=[];
if nargin>=2
    folder=varargin{1};
end
if isempty(folder)
    folder=fullocmatfile(getocmatfolder(odeObj,'userdata'));
end
fidwrite=remark(modelname(odeObj),folder,basicextension(odeObj),varargin{2:end});
if nargout==1
    varargout{1}=fidwrite;
end