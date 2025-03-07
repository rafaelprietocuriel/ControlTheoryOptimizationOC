function zeta = adj_fq(x,t,xi,theta,uy,vzq,varargin)
persistent u y z v q nx ny nz nv nq aazq aaxq aaxyq dx
  % get control w
  w = varargin{1};
  % get options
  opt = varargin{2};
  % get dimentions
  nx = length(x);
  nu = opt.Nu;
  ny = opt.Ny;
  %np = opt.Np;
  nz = opt.Nz;
  nv = opt.Nv;
  nq = opt.Nq;
  %nw = opt.Nw;
  % get parameters
  u = uy(:,1:nu);
  y = uy(:,nu+1:nu+ny);
  %p = uyp(:,nu+ny+1:nu+ny+np);
  v = vzq(1:nv);
  z = vzq(nv+1:nv+nz);
  q = vzq(nv+nz+1:nv+nz+nq);
  % get functions
  zeta = zeros(nq,1); % start constructing zeta
if isfield(opt,'funZeta')
   zeta = feval(opt.funZeta,x,t,y,z,q,xi,theta,u,v,w,varargin{2:end});
else % fzeta
  dx = diff(x);
  if isfield(opt,'funS_q')
   if ~isempty(opt.funS_q) 
    aazq = zeros(nz,nq); % auxiliary array of Jacobian matrix
    if isnumeric(opt.funS_q), aazq(:,:) = opt.funS_q;
    else aazq(:,:) = feval(opt.funS_q,x,t,y,z,v,varargin{2:end});
    end
    zeta = zeta + aazq' * theta; % Jacobian matrix is transposed !
   end
  end
  % integral part
  if isfield(opt,'funL_q')
   if ~isempty(opt.funL_q) 
    if isnumeric(opt.funL_q), zeta = zeta + trapzint(dx,opt.funL_q);
    else zeta = zeta + trapzint(dx,feval(opt.funL_q,x,t,y,z,v,w,varargin{2:end}));
    end
   end
  end
  if isfield(opt,'funF_q')
   if ~isempty(opt.funF_q) % size(funF_q) = [nx,ny,nq]
    aaxyq = zeros(nx,ny,nq); % auxiliary array
    if isnumeric(opt.funF_q), aaxyq(:,:,:) = opt.funF_q;
    else aaxyq(:,:,:) = feval(opt.funF_q,x,t,y,z,q,u,v,varargin{2:end});
    end
    aaxq = zeros(nx,nq); % auxiliary array
    for i=1:nq % for all y takes vector sum along ny
      aaxq(:,i) = aaxq(:,i) + sum(xi .* aaxyq(:,:,i),2);
    end
    zeta = zeta + trapzint(dx,aaxq);
   end
  end
end
