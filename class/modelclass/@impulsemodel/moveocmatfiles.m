function varargout=moveocmatfiles(ocObj,modelfiles,varargin)

varargout=cell(1,nargout);
[varargout{:}]=moveocmatfiles(loadmodeldata(modelname(ocObj)),modelfiles,varargin{:});