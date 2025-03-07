function out=distance(ppdePrim1,ppdePrim2,varargin)

out=norm(dependentvar(ppdePrim1)-dependentvar(ppdePrim2),varargin{:});