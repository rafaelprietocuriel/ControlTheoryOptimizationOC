function modelfiles=makefile4ocmat(ocObj,varargin)

load(ocObj,'modeldata');
modelfiles=makefile4ocmat(ocStruct,varargin{:});