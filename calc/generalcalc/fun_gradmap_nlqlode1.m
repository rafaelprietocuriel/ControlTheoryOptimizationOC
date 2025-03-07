function [u,v,w,Val,du,dv,dw,y,z,q,xi,theta,zeta] = fun_gradmap_nlqlode1(u,v,w,gmopt,varargin)
persistent fq hx ht aaqvt aav
Val = NaN;
% get options
opt = varargin{1};
% get initial conditions for primal system y
if isfield(opt,'funY0')
  if ~isempty(opt.funY0)
    if isnumeric(opt.funY0), ic = opt.funY0;
    else ic = feval(opt.funY0,w,varargin{:}); 
    end
  else error('vector funY0 must be nonempty')
  end
else error('function or vector funY0 must be given')
end
% get initial conditions for primal system z
if isfield(opt,'funZ0')
  if ~isempty(opt.funZ0)
    if isnumeric(opt.funZ0), icz = opt.funZ0;
    else icz = feval(opt.funZ0,w,varargin{:}); 
    end
  else error('vector funZ0 must be nonempty')
  end
else error('function or vector funZ0 must be given')
end
%x = opt.xmesh;
%t = opt.tspan;
% get function for nonlocal variables.
if isfield(opt,'funQ'), fq = opt.funQ;
else fq = @primal_fq_nlqlode; % if there is no any it uses default function 
end                           % that calls opt.funH
% solve primal system
[y,z,q] = nlqlode1(1,opt.funF,ic,icz,opt.xmesh,opt.tspan,fq,u,v,varargin{:});
%gnlqlode1(d,fode,ic,icz,xmesh,t,fq,u,v,varargin)
if isfield(opt,'MixConstr') % mixed constraints to control depending on state variables
    [u,v,w] = feval(opt.MixConstr,opt.xmesh,opt.tspan,y,z,q,u,v,w,varargin{:});
end
% get initial conditions for adjoint system
if isfield(opt,'funl_y')
  if ~isempty(opt.funl_y)
    if isnumeric(opt.funl_y), ic = opt.funl_y;
    else ic = feval(opt.funl_y,y(:,:,end),varargin{:}); 
    end
  else error('vector funl_y must be nonempty')
  end
else error('function or vector funl_y must be given')
end
if isfield(opt,'funl_z')
  if ~isempty(opt.funl_z)
    if isnumeric(opt.funl_z), icz = opt.funl_z;
    else icz = feval(opt.funl_z,w,varargin{:}); 
    end
  else error('vector funl_z must be nonempty')
  end
else error('function or vector funl_z must be given')
end
% calculate value function
if isfield(gmopt,'CalcVal')
  if gmopt.CalcVal
   if isfield(opt,'funJ')
    Val = feval(opt.funJ,opt.xmesh,opt.tspan,y,z,q,u,v,w,varargin{:});
   else
    hx = diff(x);
    ht = diff(t);
    Val = 0;
    if isfield(opt,'funl')
      if ~isempty(opt.funl)
        Val = Val + trapzint1q(hx,feval(opt.funl,y(:,:,end),varargin{:})); 
      end
    end
    if isfield(opt,'funLa')
      if ~isempty(opt.funLa)
        Val = Val + trapzint1q(hx,feval(opt.funLa,opt.xmesh,w,varargin{:})); 
      end
    end
    if isfield(opt,'funLt')
      if ~isempty(opt.funLt)
        Val = Val + trapzint1q(ht,feval(opt.funLt,opt.tspan,q,v,varargin{:})); 
      end
    end
    if isfield(opt,'funLta')
      if ~isempty(opt.funLta)
        Val = Val + trapzint1q(ht,trapzint2q(hx,feval(opt.funLta,opt.xmesh,opt.tspan,y,z,q,u,v,w,varargin{:}))); 
      end
    end
   end
%     if isfield(gmopt,'DisplayVal')
%       if gmopt.DisplayVal, Val
%       end
%     end
  end
end
% calculate gradient of the hamiltonian
if isfield(gmopt,'CalcGrad')
  if gmopt.CalcGrad == 0 
      %du = Val; 
      return
  end
% by default calculate gradient % else return
end
% solve adjoint system regarding concatenated u,y along dimention 2 
% and concatenated v,z,q along the first dimention. control w is
% passed as an additional parameter
[xi,theta,zeta] = nlqlode1(-1,@adj_fnlode,ic,icz,opt.xmesh,opt.tspan,@adj_fq,cat(2,u,y),cat(1,v,z,q),w,varargin{:});
% calculates the gradient of Hamiltonians
if isfield(opt,'funGradHam') % calculates the gradient by user defined function
  [du,dv,dw] = feval(opt.funGradHam,opt.xmesh,opt.tspan,y,z,q,xi,theta,zeta,u,v,w,varargin{:});
else % constract the gradiant of Hamiltonians by derivatives of the functions
  hx = diff(x);
  nv = opt.Nv;
  nz = opt.Nz;
  nt = opt.Nt;
 if ~isfield(opt,'funL_u') && ~isfield(opt,'funF_u') && ~isfield(opt,'funH_u'), du = [];
 else
  du = zeros(nx,nu,nt);
  if isfield(opt,'funL_u')    
    if isnumeric(opt.funL_u) % size(funL_u) = [nx,nu] 
      if ~isempty(opt.funL_u), error('funL_u must be nonempty') 
      end
      if isscalar(opt.funL_u), du(:) = opt.funL_u;
      elseif ndims(opt.funL_u) == 2 
        for it = 1:nt 
          du(:,:,it) = opt.funL_u; 
        end
      else du(:,:,:) = opt.funL_u; % ndims(opt.funL_u) == 3
      end
    else
      for it = 1:nt
        du(:,:,it) = feval(opt.funL_u,opt.xmesh,opt.tspan(it), ... 
          y(:,:,it),z(:,it),q(:,it),u(:,:,it),v(:,it),w,varargin{:});
      end
    end   
  end
  if isfield(opt,'funF_u') % size(funF_u) = [nx,ny,nu]
    if isnumeric(opt.funF_u) 
      if ~isempty(opt.funF_u), error('funF_u must be nonempty')
      end
      du = du + array3prod2(opt.funF_u,xi);
    else 
      for it = 1:nt
        du(:,:,it) = du(:,:,it) + array3prod2(feval(opt.funF_u, ... 
          opt.xmesh,opt.tspan(it),y(:,:,it),z(:,it),q(:,it),u(:,:,it),v(:,it),...
          varargin{:}),xi);
      end
    end
  end
  if isfield(opt,'funH_u')% size(funH_y) = [nx,nq,nu]
    if isnumeric(opt.funH_u) 
      if ~isempty(opt.funH_u) , error('funH_u must be nonempty')
      end
      if ndims(opt.funH_u) == 3 
        for iu = 1:nu 
          du(:,iu,:) = du(:,iu,:) + opt.funH_u(:,:,iu) * zeta; 
        end
      else % ndims(opt.funH_u) == 4 i.e size(opt.funH_u) == [nx,nq,nu,nt]
        for iu = 1:nu 
          for it = 1:nt 
            du(:,iu,it) = du(:,iu,it) + opt.funH_u(:,:,iu,it) * zeta(:,it); 
          end
        end
      end
    else
      aaxqu = zeros(nx,nq,nu); % auxiliary array
      for it = 1:nt 
        aaxqu(:,:,:) = feval(opt.funH_u,opt.xmesh,opt.tspan(it), ...
          y(:,:,it),z(:,it),q(:,it),u(:,:,it),v(:,it),varargin{:});
        for iu = 1:nu
          du(:,iu,it) = du(:,iu,it) + aaxqu(:,:,iu) * zeta(:,it); 
        end
      end
    end
  end
 end
 if ~isfield(opt,'funS_v') && ~isfield(opt,'funH_v') && ~isfield(opt,'funL_v') && ~isfield(opt,'funF_v'), dv = [];
 else
  dv = zeros(nv,nt); % 
  if isfield(opt,'funS_v')
    if isnumeric(opt.funS_v)
      if ~isempty(opt.funS_v), error('funS_v must be nonempty')
      end
      if ndims(opt.funS_v) == 2, dv(:,:) = opt.funS_v' * theta;
      else % if ndims(opt.funS_v) == 3 i.e. size(opt.funS_v) = [nz,nv,nt]
        for iv = 1:nv 
          dv(iv,:) = sum(theta .* reshape(opt.funS_v(:,iv,:),[nz,nt]),1); 
        end
      end
    else 
      aazv = zeros(nz,nv); % auxiliary array of Jacobian matrix
      for it = 1:nt 
        aazv(:,:) = feval(opt.funS_v,opt.tspan(it),z(:,it),q(:,it),v(:,it),varargin{:});
        dv(:,it) = dv(:,it) + aazv' * theta(:,it); % Jacobian matrix is transposed !
      end
    end    
  end
  if isfield(opt,'funH_v') % zeta * int_0_omega h_v dx
    if isnumeric(opt.funH_v) % size(funH_v) = [nx,nq,nv]
      if ~isempty(opt.funH_v), error('funH_v must be nonempty') 
      end
      if ndims(opt.funH_v) == 3, dv = dv + trapzint(hx,opt.funH_v)' * zeta; 
      else % ndims(opt.funH_v) == 4 i.e size(opt.funH_v) == [nx,nq,nv,nt]
        aaqvt = trapzint(hx,opt.funH_v); % auxiliary array
        for iv = 1:nv 
          dv(iv,:)=dv(iv,:)+sum(reshape(aaqvt(:,iv,:),[nq,nt]) .* zeta,1); 
        end
      end
    else
      for it = 1:nt 
        dv(:,it) = dv(:,it) + trapzint(hx,feval(opt.funH_v, ...
          opt.xmesh,opt.tspan(it),y(:,:,it),z(:,it),q(:,it),u(:,:,it), ...
          v(:,it),varargin{:}))' * zeta(:,it); 
      end
    end
  end
  if isfield(opt,'funL_v') % integral int_0^omega Lv dx
    if isnumeric(opt.funL_v)
      if ~isempty(opt.funL_v), error('funL_v must be nonempty') 
      end
      if size(opt.funL_v,3) == nt, dv = dv + trapzint(hx,opt.funL_v);
      elseif size(opt.funL_v,1) == nx  % if size(opt.funL_v) = [nx,nv,1]
        aav = trapzint(hx,opt.funL_v); % auxiliary array
        for iv = 1:nv 
          dv(iv,:) = dv(iv,:) +  aav(iv); 
        end
      else % if size(opt.funL_v) = [nv,1]
        for iv = 1:nv 
          dv(iv,:) = dv(iv,:) +  opt.funL_v(iv); 
        end
      end
    else 
      for it = 1:nt 
        dv(:,it) = dv(:,it) + trapzint(hx,feval(opt.funL_v, ...
          opt.xmesh,opt.tspan(it),y(:,:,it),z(:,it),q(:,it), ...
          v(:,it),u(:,:,it),w,varargin{:}));
      end
    end
  end
  if isfield(opt,'funF_v')% integral int_0^omega xi*F_v dx
    if isnumeric(opt.funF_v)% size(funF_z) = [nx,ny,nv]  
      if ~isempty(opt.funF_v), error('funF_v must be nonempty') 
      end
      dv = dv + trapzint(hx,array3prod2(opt.funF_v,xi));
    else
      for it = 1:nt 
        dv(:,it) = dv(:,it) + trapzint(hx,array3prod2(feval(opt.funF_v, ...
          opt.xmesh,opt.tspan(it),y(:,:,it),z(:,it),q(:,it),u(:,:,it), ...
          v(:,it),varargin{:}),xi(:,:,it))); 
      end
    end
  end
 end
 if ~isfield(opt,'funL_w') && ~isfield(opt,'funY0_w'), dw = [];
 else
  ht = diff(t);
  dw = zeros(nx,nw); % 
  if isfield(opt,'funL_w') % integral int_0^T L_w dt
    if isnumeric(opt.funL_w)
      if ~isempty(opt.funL_w), error('funL_w must be nonempty') 
      end
      if size(opt.funL_w,3) == nt, dw = dw + trapzint(ht,opt.funL_v,3);
      elseif size(opt.funL_w,1) == nx  % if size(opt.funL_w) = [nx,nw,1]
        dw = dw + opt.funL_w;
      else % if size(opt.funL_w) = [nw,1]
        for iw = 1:nw 
          dw(:,iw) = dw(:,iw) +  opt.funL_w(iw); 
        end
      end
    else % size(opt.funL_w) = [nx,nw,1]
      for ix = 1:nx 
        dw(ix,:) = dw(ix,:) + trapzint(ht,feval(opt.funL_w, ...
          opt.xmesh,opt.tspan(it),reshape(y(ix,:,:),[ny,nt]),z,q, ...
          v,reshape(u(ix,:,:),[nu,nt]),w(ix,:),varargin{:}),3);
      end
    end
  end
  if isfield(opt,'funY0_w')% scalar multiplication int_0^omega xi*Y0_w dx
    if isnumeric(opt.funY0_w)% size(funY0_w) = [nx,ny,nw]  
      if ~isempty(opt.funY0_w), error('funY0_w must be nonempty') 
      end
      dw = dw + array3prod2(opt.funY0_w,xi(:,:,1));
    else
      for ix = 1:nx 
        dw = dw + array3prod2(feval(opt.funY0_w, ...
          opt.xmesh,w,varargin{:}),xi(:,:,1));
      end
    end
   end
  end
 end
end