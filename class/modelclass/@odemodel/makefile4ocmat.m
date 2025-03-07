function modelfiles=makefile4ocmat(odeObj,varargin)

load(odeObj,'modeldata');
modelfiles=makefile4ocmat(ocStruct,varargin{:});