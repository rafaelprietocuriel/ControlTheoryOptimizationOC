function [f,s] = godemd_prm(x,t,y,z,u,v,varargin)
% general non local ODE on moving domain primal function
persistent b
  % get control w
  % w = varargin{1};
  % get options
  %opt = varargin{2};
  % get functions
if isfield(opt,'funF')
  b = x <= z(1); 
  [f,s] = feval(opt.funF,x(b),t,y(b,:),z,u(b,:),v,varargin{2:end});
else error('funF must be specified!');
end