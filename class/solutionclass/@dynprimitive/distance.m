function out=distance(dynPrim1,dynPrim2,varargin)

out=norm(dependentvar(dynPrim1)-dependentvar(dynPrim2),varargin{:});