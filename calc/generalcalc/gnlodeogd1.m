function [y,z,p,q] = gnlodeogd1(d,fode,ic,icz,fbc,xmesh,t,fpq,u,v,varargin)
%   [Y,P,Q] = GNLODEOGD1(d,ODEFUN,IC,BCFUN,XMESH,TSPAN,PQFUN) solves initial-boundary
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
%   Y = BCFUN(T,Z,Q) is a function that evaluates the
%   components of the boundary conditions at scalar time T 
%   depending on vector papameter Q,  size(Q) = [nq].
%   Y is column vector, evaluated at a boundary. When d=1 it is the right boundary, 
%   while when d=-1 it is the left boundary.
%   If state variables in dimention I shouldn't have any boundary conditions 
%   then BCFUN should return Y(I) = NaN
%   If boundary conditions are constants then instead of functiona henler 
%   argument BCFUN can be set as real array of length ny !
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
%   [P,Q] = PQFUN(XMESH,T,Y) is a function that evaluates 
%   the nonlocal variables p and q. The input arguments are
%   XMESH is the vector. size(XMESH) = [nx,1], 
%   T is the scalat time
%   Y is the matrix of state variables as XMESH, size(Y) = [nx,ny]
%   PQFUN returns [P,Q], size(P)  = [nx,np], size(Q)  = [nq,1].
%
%   GNLQLPDEWODE1 returns matrices
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
%   Revision date: 2011/02/11
persistent f
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
  %x0 = nx; % index at which boundary conditions are set
  dt0 = nt - 1; % the first index for array dt of time intervals
  dtf = 1; % the last index for array dt of time intervals
  i0 = nt; % index for array t of time instances, beginning of interval dt(i) 
  nx0 = nx; %nx0 = nx;
  nxM = nx;
  nx = nx - nt; %nx = nx - nt + 1;
  zxmesh = xmesh;
case 1
  %x0 = 1; % index at which boundary conditions are set
  dt0 = 1; % the first index for array dt of time intervals
  dtf = nt - 1; % the last index for array dt of time intervals
  i0 = 1; % index for array t of time instances, beginning of interval dt(i) 
  nx0 = nx + 1; % nx0 = nx; 
  nxM = nx + nt; %nxM = nx+nt-1;
  %zxmesh = [xmesh;zeros(nt-1,size(xmesh,2))];
  zxmesh = [xmesh;zeros(nt,size(xmesh,2))];
  zxmesh(nx0,:) = icz(:)';
otherwise
  error('d must be 1 or -1.')
end
% geting the dimentions of vectors and their initial values
if isempty(ic)
  ny = 0;
  y = zeros(nxM,ny,nt); % declare the array of state variables
else
  ny = size(ic,2); % get number of state variables from IC second dimention
  y = zeros(nxM,ny,nt); % declare the array of state variables
  nc = size(ic,1);
  y(1:nc,:,i0) = ic(1:nc,:); % write the initial conditions for all x
end
% upper bound of the domain that is why size(z,1) = ny
%z = zeros(ny,nt); % declare the array of state variables
%z(:,i0) = xmesh(end,:)'; % write the initial conditions of z
nz = size(icz,1); % get number of state variables from IC first dimention
z = zeros(nz,nt); % declare the array of state variables
z(:,i0) = icz; % write the initial conditions of z

if nargin > 8 % exist(u)
  if isempty(u), u = zeros(nxM,0,nt); end
  %if size(u,1) ~= nx+nt-1, error('size(u,1) = nx+nt-1'), end
else
  u = zeros(nxM,0,nt);
end

if nargin > 9 % exist(v)
  if isempty(v), v = zeros(0,nt); end
else
  v = zeros(0,nt); 
end

if nargin > 7 %exist('fpq')
  [p0,q0] = feval(fpq,zxmesh(1:nx0,:),t(i0),y(1:nx0,:,i0),z(:,i0),u(1:nx0,:,i0),v(:,i0),varargin{:});
  if isempty(p0) 
    np = 0;
    p0 = zeros(nxM,0);
  else
    np = size(p0,2);
    %if size(p0,1) ~= nx+nt-1, error('size(p0,1) = nx+nt-1'), end
  end
  if isempty(q0), nq = 0; 
  else nq = size(q0,1);
  end
end
p = zeros(nxM,np,nt); % may be declare this variable as persistent ?????
q = zeros(nq,nt);  % may be declare this variable as persistent ?????
p(1:nx0,:,i0) = p0(1:nx0,:);
if nq ~= 0, q(:,i0) = q0;end

% arrays declaration before the loop for faster perfomance
%xymesh = xmesh * ones(1,ny); % prepare rectangular mesh for state variables
%y_int = zeros(nx,ny); % local array for integration (on rectangular mesh)
y_t = zeros(nxM,ny);% local array of values on nonrectangular mesh
z_t = zeros(nz);% local array of z(t)
fz = z_t;% slopes of domain upper bound and other variables
%xymesh_t = zeros(nx,ny);% local array of nonrectangular mesh
% local array of boundary conditions
bc_t = zeros(ny);
if isnumeric(fbc) 
  if isempty(fbc) % if there is no boundary conditions
    bc_t(:) = NaN; 
  else
    if length(fbc) ~= ny, error('length(fbc) = ny'), end
    bc_t = fbc;
  end
else % get the boundary conditions for t(i0) with q for t(i0) !
    bc_t = feval(fbc,t(i0),z(:,i0),q(:,i0),v(:,i0),varargin{:});% and first order approx. of z
end
if d == 1 % set boundary conditions at the initial moment of time
%     if any(isnan(bc_t)) == 1
%         y(nx0,~isnan(bc_t),1) = bc_t(~isnan(bc_t));
%     end
    booo = ~isnan(bc_t);
    if any(booo), y(nx0,booo,1) = bc_t(booo); end
end
%f = zeros(nx+nt,ny);% right hand side values
%ff = zeros(nx+nt,ny);% right hand side values
% the main loop over time. nt-1 times
for i = dt0:d:dtf % index for array dt of time increments
  if d == 1, dti = dt(i); else dti = -dt(i); end % set sign to the time interval
  i1 = i0 + d; % index for array t of time instances, ending of interval dt(i)
  nx1 = nx0 + d;
  [f(1:nx0,:),fz(:)] = feval(fode,zxmesh(1:nx0,:),t(i0),y(1:nx0,:,i0),z(:,i0),p(1:nx0,:,i0),q(:,i0),u(1:nx0,:,i0),v(:,i0),varargin{:}); % get functions f and c
  y_t(1:nx0,:) = y(1:nx0,:,i0) + dti * f(1:nx0,:); % get approximations of y which generally are not on the grid nodes
  z_t = z(:,i0) + dti * fz; % get approximations of z
  %y_t(isinf(y_t(1:nx0,:)) | ~isreal(y_t(1:nx0,:))) = 0;
  %z_t(isinf(z_t) | ~isreal(z_t)) = 0;
  if d == 1 
      %y_t = max(0,y_t);
      zxmesh(nx1,1:nx2) = z_t(1:nx2)';
      if ~isnumeric(fbc) % get the boundary conditions for t(i1) but with q for t(i0) !
        bc_t = feval(fbc,t(i1),z_t,q(:,i0),v(:,i1),varargin{:});% and first order approx. of z
      else bc_t = fbc;
      end
      bo = isinf(bc_t) | z(1:nx2,i0) == z_t(1:nx2);%bo = isinf(bc_t);
      bc_t(bo) = y_t(nx0,bo); % if there is a singularity use the previous nx0 point
      y_t(nx1,:) = bc_t;
  end
  % p and q
  if np ~= 0 || nq ~= 0 
    [p0,q0] = feval(fpq,zxmesh(1:nx1,:),t(i1),y_t(1:nx1,:),z_t,u(1:nx1,:,i1),v(:,i1),varargin{:}); 
    if np ~= 0, p(1:nx1,:,i1) = p0(1:nx1,:);end
    if nq ~= 0, q(:,i1) = q0;end
  end
  % get functions f and c on the end of interval dt(i)
  [f(1:nx1,:),fz(:)] = feval(fode,zxmesh(1:nx1,:),t(i1),y_t(1:nx1,:),z_t,p(1:nx1,:,i1),q(:,i1),u(1:nx1,:,i1),v(:,i1),varargin{:});
  y(1:nx0,:,i1) = 0.5 * (y_t(1:nx0,:) + y(1:nx0,:,i0) + dti * f(1:nx0,:)); % take average between two approximations of y
  z(:,i1) = 0.5 * (z_t + z(:,i0) + dti * fz); % take average between two approximations of z
  %y(isinf(y(1:nx0,:,i1)) | ~isreal(y(1:nx0,:,i1)),i1) = 0;
  %z(isinf(z(:,i0)) | ~isreal(z(:,i0)),i0) = 0;
  %z(isinf(z(:,i1)) | ~isreal(z(:,i1)),i1) = 0;
  if d == 1
    %y(1:nx0,:,i1) = max(0,y(1:nx0,:,i1));
    zxmesh(nx1,1:nx2) = z(1:nx2,i1)';
    if ~isnumeric(fbc) % get the boundary conditions for t(i1) and with q for t(i1) !
      bc_t = feval(fbc,t(i1),z(:,i1),q(:,i1),v(:,i1),varargin{:});
    else bc_t = fbc;
    end
    bo = isinf(bc_t) | z(1:nx2,i0) == z(1:nx2,i1);
    bc_t(bo) = y(nx0,bo,i1); % if there is a singularity use the previous point nx0 and time i1
    y(nx1,:,i1) = bc_t;
  end
  % nonlocal variables p and q
  if np ~= 0 || nq ~= 0 
    [p0,q0] = feval(fpq,zxmesh(1:nx1,:),t(i1),y(1:nx1,:,i1),z(:,i1),u(1:nx1,:,i1),v(:,i1),varargin{:}); 
    if np ~= 0, p(1:nx1,:,i1) = p0(1:nx1,:);end
    if nq ~= 0, q(:,i1) = q0;end
  end
  % get the boundary condition again with improved q 
  % it does not influence the order of convergance but increases accuracy
  % if boundary conditions depend on q
  if d == 1
    if ~isnumeric(fbc)
      bc_t = feval(fbc,t(i1),z(:,i1),q(:,i1),v(:,i1),varargin{:});
    else bc_t = fbc;
    end
    bo = isinf(bc_t) | z(1:nx2,i0) == z(1:nx2,i1);% | fz == 0;
    bc_t(bo) = y(nx0,bo,i1); % if there is a singularity use the previous point nx0 and time i0
    y(nx1,:,i1) = bc_t;
  end
  %
  i0 = i0 + d;     % index for array t of time instances, beginning of interval dt(i) 
  nx0 = nx0 + d;
end