function varargout=remark(dgObj,varargin)
%
% REMARK creates or opens a text file for remarks.
%
% REMARK(OCOBJ) creates a file with the basic filename of the model OCOBJ,
% see ocmat/tools/remark for further information.

folder=[];
if nargin>=2
    folder=varargin{1};
end
if isempty(folder)
    folder=fullocmatfile(getocmatfolder(dgObj,'userdata'));
end
fidwrite=remark(modelname(dgObj),folder,basicextension(dgObj));
if nargout==1
    varargout{1}=fidwrite;
end