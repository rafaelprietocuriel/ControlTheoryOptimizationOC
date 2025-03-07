function varargout=remark(ocObj,varargin)
%
% REMARK creates or opens a text file for remarks.
%
% REMARK(OCOBJ) creates a file with the basic filename of the model OCOBJ,
% see ocmat/tools/remark for further information.
%
% REMARK(OCOBJ,FOLDER) the file is saved inside the folder FOLDER. If
% FOLDER is empty the user model data folder is used.
%
% REMARK(OCOBJ,FOLDER,WD) if WD is TRUE the actual date is appended at the
% end of the file, otherwise not.

folder=[];
if nargin>=2
    folder=varargin{1};
end
if isempty(folder)
    folder=fullocmatfile(getocmatfolder(ocObj,'userdata'));
end
fidwrite=remark(modelname(ocObj),folder,basicextension(ocObj),varargin{2:end});
if nargout==1
    varargout{1}=fidwrite;
end