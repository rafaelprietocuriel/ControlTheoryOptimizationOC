function f = adj_fode(t,theta,vz,varargin)
persistent v z nz aazz
  % get control z_0
  z_0 = varargin{1};
  % get options
  opt = varargin{2};
  % get dimentions
  %nu = opt.Nu;
  nv = opt.Nv;
  nz = opt.Nz;
  %nv = opt.Nv;
  % get parameters
  v = vz(1:nv);
  z = vz(nv+1:nv+nz);
  % get functions
 if isfield(opt,'fun_adj')
  f = feval(opt.fun_adj,t,z,theta,v,z_0,varargin{2:end});
  return
 elseif isfield(opt,'funS_adj')
  f = feval(opt.funS_adj,t,z,theta,v,z_0,varargin{2:end});
 else
  f = zeros(nz,1);
  if isfield(opt,'funS_z')
   if ~isempty(opt.funS_z) 
    aazz = zeros(nz,nz); % auxiliary array of Jacobian matrix
    if isnumeric(opt.funS_z), aazz(:,:) = opt.funS_z;
    else aazz(:,:) = feval(opt.funS_z,t,z,v,varargin{2:end});
    end
    f = f + aazz' * theta; % Jacobian matrix is transposed !
   end
  end
  if isfield(opt,'funL_z')
   if ~isempty(opt.funL_z) 
    if isnumeric(opt.funL_z), f = f + opt.funL_z;
    else f = f + feval(opt.funL_z,t,z,v,z_0,varargin{2:end});
    end
   end
  end
 end
 f = - f;

  