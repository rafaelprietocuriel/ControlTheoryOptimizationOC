function [fxi,c,ftheta] = adj_fqnlpdewode(x,t,xi,theta,eta,zeta,uyp,vzq,varargin)
persistent u y z v p q nx ny nz nv np nq aaxyy aaxqy aaxyz aaxz dx
  % get control w
  w = varargin{1};
  % get options
  opt = varargin{2};
  % get dimentions
  nx = length(x);
  nu = opt.Nu;
  ny = opt.Ny;
  np = opt.Np;
  nz = opt.Nz;
  nv = opt.Nv;
  nq = opt.Nq;
  %nw = opt.Nw;
  % get parameters
  u = uyp(:,1:nu);
  y = uyp(:,nu+1:nu+ny);
  p = uyp(:,nu+ny+1:nu+ny+np);
  v = vzq(1:nv);
  z = vzq(nv+1:nv+nz);
  q = vzq(nv+nz+1:nv+nz+nq);
  % get functions
if isfield(opt,'fun_adj')
  [fxi,c,ftheta] = feval(opt.fun_adj,x,t,y,z,p,q,xi,theta,eta,zeta,u,v,w,varargin{2:end});
else
 % fxi
 dx = diff(opt.xmesh);
 if isfield(opt,'funF_adj')
  fxi = feval(opt.funF_adj,x,t,y,z,p,q,xi,theta,eta,zeta,u,v,w,varargin{2:end});
 else
  fxi = zeros(nx,ny);
  if isfield(opt,'funL_y')
   if ~isempty(opt.funL_y)
    if isnumeric(opt.funL_y), fxi = fxi + opt.funL_y;
    else fxi = fxi + feval(opt.funL_y,x,t,y,z,p,q,u,v,w,varargin{2:end});
    end
   end
  end
  if isfield(opt,'funF_y')
   if ~isempty(opt.funF_y) % size(funF_y) = [nx,ny,ny]
    aaxyy = zeros(nx,ny,ny); % auxiliary array
    if isnumeric(opt.funF_y), aaxyy(:,:,:) = opt.funF_y;
    else aaxyy(:,:,:) = feval(opt.funF_y,x,t,y,z,p,q,u,v,varargin{2:end});
    end
    for i=1:ny % for all y takes vector sum along ny
      fxi(:,i) = fxi(:,i) + sum(xi .* aaxyy(:,:,i),2);
    end
   end
  end
  if isfield(opt,'funH_y')
   if ~isempty(opt.funH_y) % size(funH_y) = [nx,nq,ny]
    aaxqy = zeros(nx,nq,ny); % auxiliary array
    if isnumeric(opt.funH_y), aaxqy(:,:,:) = opt.funH_y;
    else aaxqy(:,:,:) = feval(opt.funH_y,x,t,y,z,p,q,u,v,varargin{2:end});
    end
    for i=1:ny % for all y takes vector sum along nq
      fxi(:,i) = fxi(:,i) + aaxqy(:,:,i) * zeta;
    end
   end
  end
 end 
 fxi = - fxi;
 % ftheta
 if isfield(opt,'funS_adj')
  ftheta = feval(opt.funS_adj,x,t,y,z,p,q,theta,v,w,varargin{2:end});
 else
  ftheta = zeros(nz,1);
  if isfield(opt,'funS_z')
   if ~isempty(opt.funS_z) 
    aazz = zeros(nz,nz); % auxiliary array of Jacobian matrix
    if isnumeric(opt.funS_z), aazz(:,:) = opt.funS_z;
    else aazz(:,:) = feval(opt.funS_z,t,z,q,v,varargin{2:end});
    end
    ftheta = ftheta + aazz' * theta; % Jacobian matrix is transposed !
   end
  end
  % integral part
  if isfield(opt,'funL_z')
   if ~isempty(opt.funL_z) 
    if isnumeric(opt.funL_z), ftheta = ftheta + trapzint(dx,opt.funL_z);% L_z could be integrated in advance !
    else ftheta = ftheta + trapzint(dx,feval(opt.funL_z,x,t,y,z,p,q,v,w,varargin{2:end}));
    end
   end
  end
  if isfield(opt,'funF_z')
   if ~isempty(opt.funF_z) % size(funF_z) = [nx,ny,nz]
    aaxyz = zeros(nx,ny,nz); % auxiliary array
    if isnumeric(opt.funF_z), aaxyz(:,:,:) = opt.funF_z;
    else aaxyz(:,:,:) = feval(opt.funF_z,x,t,y,z,p,q,u,v,varargin{2:end});
    end
    aaxz = zeros(nx,nz); % auxiliary array
    for i=1:nz % for all y takes vector sum along ny
      aaxz(:,i) = aaxz(:,i) + sum(xi .* aaxyz(:,:,i),2);
    end
    ftheta = ftheta + trapzint(dx,aaxz);
   end
  end
 end
 ftheta = - ftheta;
end

  