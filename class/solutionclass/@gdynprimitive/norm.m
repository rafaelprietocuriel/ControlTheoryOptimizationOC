function out=norm(gdynPrim,varargin)

X=dependentvar(gdynPrim);
out=norm(X{1},varargin{:});