function varargout=moveocmatfiles(dgObj,modelfiles,varargin)

varargout=cell(1,nargout);
if ~isempty(modelfiles)
    testfile=fullocmatfile(getocmatfolder('out'),[modelname(dgObj) modelfiles(1).name '.m']);
    if ~strcmp(testfile,modelfiles(1).filename)
        ocmatmsg('Model ''%s'' and files are not consistent.\n',modelname(dgObj))
        return
    end
end
[varargout{:}]=moveocmatfiles(loadmodeldata(modelname(dgObj)),modelfiles,varargin{:});