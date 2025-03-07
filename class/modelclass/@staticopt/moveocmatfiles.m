function varargout=moveocmatfiles(ocObj,modelfiles,varargin)

varargout=cell(1,nargout);
if ~isempty(modelfiles)
    testfile=fullocmatfile(getocmatfolder('out'),[modelname(ocObj) modelfiles(1).name '.m']);
    if ~strcmp(testfile,modelfiles(1).filename)
        ocmatmsg('Model ''%s'' and files are not consistent.\n',modelname(ocObj))
        return
    end
end
[varargout{:}]=moveocmatfiles(loadmodeldata(modelname(ocObj)),modelfiles,varargin{:});