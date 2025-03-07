function modelfiles=makefile4ocmat(dgObj,varargin)

load(dgObj,'modeldata');
modelfiles=makefile4ocmat(ocStruct,varargin{:});