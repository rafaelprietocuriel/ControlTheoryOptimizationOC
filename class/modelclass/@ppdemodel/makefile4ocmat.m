function modelfiles=makefile4ocmat(ppdeObj,varargin)

load(ppdeObj,'modeldata');
modelfiles=makefile4ocmat(ocStruct,varargin{:});