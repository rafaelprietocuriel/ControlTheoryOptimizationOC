function [y,z,p,q] = gnlode1(d,fode,ic,icz,xmesh,t,fpqy0,u,v,varargin)
%   [Y,P,Q] = GNLODE1(d,ODEFUN,IC,BCFUN,XMESH,TSPAN,QPFUN) solves initial-boundary
%   value problems for small systems of nonlocal quasi-linear ODEs in one
%   space variable x and time t along with system of ODEs.    
%
%   Dy/Dt = f(x,t,y(x,t),p(x,t),q(t))
%
%   Dz(t)/Dt = fz(x,t,y(x,t),z(t),q(t),p(x,t)) 
%
%   where vector f(x,t,y(x,t),p(x,t),q(t)) 
%   can depend on nonlocal functions with respect to space x
%
%   p(x,t) = fp(x,t,y(x,t),p(x,t),q(t)) 
%
%   q(t) = fq(x,t,y(x,t),p(x,t),q(t),p(x,t))
%
%   functions fp and fq can contain any integral trasformations along
%   the space x
%
%   Parameter d must be -1 or 1. It sets the direction of the problem
%   (positive or negative) simultaniously both in time t and argument x. 
%   The PDEs hold for t0 <= t <= tf and x0 <= x <= xf when d = 1 
%   or t0 >= t >= tf and x0 >= x >= xf when d = -1
%
%   All multi-dimentional arrays have the following structure rule:
%   The first dimention is space variable index if any.
%   The second dimention is state variable index if any and if there is
%   the space index as the first dimention.
%   The last dimention is time index if any.
%
%   [F,C,FZ] = PDEFUN(XMESH,T,Y,Z,P,Q) is a function that evaluates the quantities
%   defining the differential equation. The input arguments are
%   XMESH is the vector of values x
%   T is the scalat time
%   Y is the matrix of state variables. size(Y) = [nx,ny]
%   P is the matrix. size(P) = [nx,np]
%   Q is the vector. size(Q) = [nq]
%   PDEFUN returns matrices
%   F is the right hand side of the PDE. size(F) = [nx,ny]
%   C contains coefficients at detivatives Dy/Dx. size(C) = [nx,ny]
%   it determines the direction of caracteristics of the equation
%   If C is empty array then it is the same as C=0, i.e. we solve ODEs
%   FZ is the right hand side of the ODE. size(FZ) = [nz,1]
%
%   ICY is the matrix containing the initial conditions. size(ICY) = [nx,ny]
%   ICZ is the vector containing the initial conditions. size(ICZ) = [nz,1]
%
%   XMESH must be strictly increasing array of values x at which the
%   solution is evaluetad
%
%   TSPAN must be strictly increasing array of time inctances at which the
%   solution is evaluetad
%
%   NLQLDE1 returns values of the solution on a mesh provided as the input
%   array XMESH. The entries of XMESH must satisfy 
%       a = XMESH(1) < XMESH(2) < ... < XMESH(NX) = b 
%   for some NX >= 3.  The entries of TSPAN must satisfy 
%       t0 = TSPAN(1) < TSPAN(2) < ... < TSPAN(NT) = tf 
%   for some NT >= 3. 
%  
%   [P,Q,Y0] = PQFUN(XMESH,T,Y) is a function that evaluates 
%   the nonlocal variables p and q. The input arguments are
%   XMESH is the vector. size(XMESH) = [nx,1], 
%   T is the scalat time
%   Y is the matrix of state variables as XMESH, size(Y) = [nx,ny]
%   PQFUN returns [P,Q,Y0], size(P)  = [nx,np], size(Q)  = [nq,1] , size(Y0).
%
%   GNLODE1 returns matrices
%   Y is matrix of state variables y, size(Y) = [nx,ny,nt]
%   Z is matrix of variable z, size(Z) = [nz,nt]
%   P is matrix of variable p, size(P) = [nx,np,nt]
%   Q is matrix of variable q, size(Q) = [nq,nt]
%
%   GNLQLPDEWODE1 uses Heun's method (also called the modified Euler's method)
%   to solve ODE along characteristics calculated on each time step.
%   Quick linear interpolation function INTERPLQ is used to obtain values for the
%   same XMESH on each time step.
%
%   Vladimir M. Veliov
%   Anton O. Belyakov
%   Revision date: 2010/11/02
persistent f fz y_t z_t
% check input parameters
if nargin < 6 
  error('Not enough input arguments')
end

%t = t(:); % make shure that t is a one-dimentional array
nt = length(t);
if nt < 3
  error('TSPAN must have at least 3 entries.')
end
dt = diff(t); % time increments
% if any(dt <= 0)
%   error('The entries of TSPAN must be strictly increasing.')
% end

%xmesh = xmesh(:); % make shure that xmesh is a one-dimentional array
% nz1 = size(xmesh,2); 
% only first nx2 elements of vector z are the bounds
[nx, nx2] = size(xmesh);
if nx < 3
  error('XMESH must have at least 3 entries.')
end
% h = diff(xmesh); % argument increments
% if any(h <= 0)
%   error('The entries of XMESH must be strictly increasing.')
% end
% define initial and final indices of XMESH, TSPAN, x, and dt
switch d % direction of the problem forward (d=1) or backward (d=-1)
case -1 
  dt0 = nt - 1; % the first index for array dt of time intervals
  dtf = 1; % the last index for array dt of time intervals
  i0 = nt; % index for array t of time instances, beginning of interval dt(i) 
case 1
  dt0 = 1; % the first index for array dt of time intervals
  dtf = nt - 1; % the last index for array dt of time intervals
  i0 = 1; % index for array t of time instances, beginning of interval dt(i) 
otherwise
  error('d must be 1 or -1.')
end
% geting the dimentions of vectors and their initial values
if isempty(ic)
  ny = 0;
  y = zeros(nx,ny,nt); % declare the array of state variables
else
  ny = size(ic,2); % get number of state variables from IC second dimention
  y = zeros(nx,ny,nt); % declare the array of state variables
  nc = size(ic,1);
  y(1:nc,:,i0) = ic(1:nc,:); % write the initial conditions for first x(1:nc)
end
nz = size(icz,1); % get number of state variables from IC first dimention
z = zeros(nz,nt); % declare the array of state variables
if nz > 0, z(:,i0) = icz; end % write the initial conditions of z

if nargin > 7 % exist(u)
  if isempty(u), u = zeros(nx,0,nt); end
else
  u = zeros(nx,0,nt);
end

if nargin > 8 % exist(v)
  if isempty(v), v = zeros(0,nt); end
else
  v = zeros(0,nt); 
end

np = 0;
nq = 0;
ny0 = 0;
if nargin > 6 %exist('fpq')
  switch nargout(fpqy0)
    case 1
      [p0] = feval(fpqy0,xmesh,t(i0),y(:,:,i0),z(:,i0),u(:,:,i0),v(:,i0),varargin{:});
    case 2
      [p0,q0] = feval(fpqy0,xmesh,t(i0),y(:,:,i0),z(:,i0),u(:,:,i0),v(:,i0),varargin{:});
      nq = size(q0,1);
    case 3
      [p0,q0,y0] = feval(fpqy0,xmesh,t(i0),y(:,:,i0),z(:,i0),u(:,:,i0),v(:,i0),varargin{:},i0,z(:,i0));
      nq = size(q0,1);
      ny0 = size(y0,2);
      % write the initial conditions for first x(1:nc)
      if ny0 ~= 0 && nc < nx 
          y(nc+1:end,:,i0) = y0(nc+1:end,:); 
      end
  end
  np = size(p0,2);
end
p = zeros(nx,np,nt); % may be declare this variable as persistent ?????
q = zeros(nq,nt);  % may be declare this variable as persistent ?????
if np ~= 0, p(:,:,i0) = p0(:,:);end
if nq ~= 0, q(:,i0) = q0;end

% arrays declaration before the loop for faster perfomance
%y_t = zeros(nx,ny);% local array of values on nonrectangular mesh
%z_t = zeros(nz);% local array of z(t)
%fz = z_t;% slopes of domain upper bound and other variables
% the main loop over time. nt-1 times
for i = dt0:d:dtf % index for array dt of time increments
  if d == 1, dti = dt(i); else dti = -dt(i); end % set sign to the time interval
  i1 = i0 + d; % index for array t of time instances, ending of interval dt(i)
  % get functions f and c
  if nz > 0 
    [f,fz] = feval(fode,xmesh,t(i0),y(:,:,i0),z(:,i0),p(:,:,i0),q(:,i0),u(:,:,i0),v(:,i0),varargin{:}); 
    z_t = z(:,i0) + dti * fz; % get approximations of z
  else
    [f] = feval(fode,xmesh,t(i0),y(:,:,i0),z(:,i0),p(:,:,i0),q(:,i0),u(:,:,i0),v(:,i0),varargin{:}); 
  end
  y_t = y(:,:,i0) + dti * f(:,:); % get approximations of y which generally are not on the grid nodes
  % p and q and y0
  if ny0 ~= 0
    [p0,q0,y0] = feval(fpqy0,xmesh(:,:),t(i1),y_t(:,:),z_t,u(:,:,i1),v(:,i1),varargin{:},i1,z(:,i0));   
    y_t(:,:) = y0;
  elseif nq ~= 0
    [p0,q0] = feval(fpqy0,xmesh(:,:),t(i1),y_t(:,:),z_t,u(:,:,i1),v(:,i1),varargin{:});
  elseif np ~= 0
    [p0] = feval(fpqy0,xmesh(:,:),t(i1),y_t(:,:),z_t,u(:,:,i1),v(:,i1),varargin{:});
  end
  if np ~= 0, p(:,:,i1) = p0(:,:);end
  if nq ~= 0, q(:,i1) = q0;end
  % get functions f and c on the end of interval dt(i)
  if nz > 0 
    [f,fz] = feval(fode,xmesh,t(i1),y_t(:,:),z_t,p(:,:,i1),q(:,i1),u(:,:,i1),v(:,i1),varargin{:});
    z(:,i1) = 0.5 * (z_t + z(:,i0) + dti * fz); % take average between two approximations of z
  else
    [f] = feval(fode,xmesh,t(i1),y_t(:,:),z_t,p(:,:,i1),q(:,i1),u(:,:,i1),v(:,i1),varargin{:});
  end
  y(:,:,i1) = 0.5 * (y_t(:,:) + y(:,:,i0) + dti * f(:,:)); % take average between two approximations of y
  % nonlocal variables p and q
  if np ~= 0 || nq ~= 0 
    [p0,q0] = feval(fpqy0,xmesh(:,:),t(i1),y(:,:,i1),z(:,i1),u(:,:,i1),v(:,i1),varargin{:}); 
    if np ~= 0, p(:,:,i1) = p0(:,:);end
    if nq ~= 0, q(:,i1) = q0;end
  end
  if ny0 ~= 0
    [p0,q0,y0] = feval(fpqy0,xmesh(:,:),t(i1),y(:,:,i1),z(:,i1),u(:,:,i1),v(:,i1),varargin{:},i1,z(:,i0));    
    y(:,:,i1) = y0;
  elseif nq ~= 0
    [p0,q0] = feval(fpqy0,xmesh(:,:),t(i1),y(:,:,i1),z(:,i1),u(:,:,i1),v(:,i1),varargin{:});
  elseif np ~= 0
    [p0] = feval(fpqy0,xmesh(:,:),t(i1),y(:,:,i1),z(:,i1),u(:,:,i1),v(:,i1),varargin{:});
  end
  if np ~= 0, p(:,:,i1) = p0(:,:);end
  if nq ~= 0, q(:,i1) = q0;end
  %
  i0 = i0 + d;     % index for array t of time instances, beginning of interval dt(i) 
end