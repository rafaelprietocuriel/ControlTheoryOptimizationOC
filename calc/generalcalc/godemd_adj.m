function [fxi,ftheta] = godemd_adj(x,t,xi,theta,uy,vz,varargin)
% general non local ODE on moving domain adjoint function
persistent u y z v ny nz nv
  % get control w
  w = varargin{1};
  % get options
  opt = varargin{2};
  % get dimentions
  nx = length(x);
  nu = opt.Nu;
  ny = opt.Ny;
  nz = opt.Nz;
  nv = opt.Nv;
  %nw = opt.Nw;
  % get parameters
  u = uy(:,1:nu);
  y = uy(:,nu+1:nu+ny);
  v = vz(1:nv);
  z = vz(nv+1:nv+nz);
  % get functions
if isfield(opt,'fun_adj')
  [fxi,ftheta] = feval(opt.fun_adj,x,t,y,z,xi,theta,u,v,w,varargin{2:end});
else error('fun_adj must be specified!');
end

  