function [p,q] = primal_fpq(x,t,y,u,v,varargin)
  % get options
  opt = varargin{1};
  % get dimentions
  nx = length(x);
  % get functions
  h = diff(x);
  if isfield(opt,'funG')
      p = trapzint(h,feval(opt.funG,x,t,y,u,v,varargin{:}));
  else p = zeros(nx,0);
  end
  if isfield(opt,'funH')
      q = trapzint(h,feval(opt.funH,x,t,y,p,u,v,varargin{:}));
  else q = [];
  end
end