function open(ocObj,filename,varargin)
%
% REMARK creates or opens a text file for remarks.
%
% REMARK(OCOBJ) creates a file with the basic filename of the model OCOBJ,
% see ocmat/tools/remark for further information.

folder=[];
if nargin>=3
    folder=varargin{1};
end
if isempty(folder)
    folder=fullocmatfile(getocmatfolder(ocObj,'specificmodel'));
end
fullfilename=fullfile(folder,[modelname(ocObj) filename '.m']);
if ~exist(fullfilename,'file')
    ocmatmsg([modelname(ocObj) filename '.m does not exist.'])
    return
end
if ispc
    try
        winopen(fullfilename);
    catch
        edit(fullfilename)
    end
else
    edit(fullfilename)
end
