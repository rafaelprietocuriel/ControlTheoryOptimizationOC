function [y,z,p,q] = gnlqlpdewlode1(d,fpde,ic,icz,fbc,xmesh,t,fpq,u,v,varargin)
%   [Y,P,Q] = GNLQLPDEWODE1(d,PDEFUN,IC,ICZ,BCFUN,XMESH,TSPAN,PQFUN) solves initial-boundary
%   value problems for small systems of nonlocal quasi-linear PDEs in one
%   space variable x and time t along with system of ODEs.    
%
%   Dy/Dt + c(x,t,y(x,t),p(x,t),q(t)) * Dy/Dx = f(x,t,y(x,t),p(x,t),q(t))
%
%   Dz(t)/Dt = fz(x,t,y(x,t),z(t),q(t),p(x,t)) 
%
%   where matrix c(x,t,p(x,t),q(t)) and vector f(x,t,y(x,t),p(x,t),q(t)) 
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
%   Revision date: 2012/08/24
persistent xb % coordinates of discontinuity which are kept to the next function call
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
nx = length(xmesh);
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
  x0 = nx; % index at which boundary conditions are set
  dt0 = nt - 1; % the first index for array dt of time intervals
  dtf = 1; % the last index for array dt of time intervals
  i0 = nt; % index for array t of time instances, beginning of interval dt(i) 
  i1 = nt - 1; % index for array t of time instances, ending of interval dt(i)
case 1
  x0 = 1; % index at which boundary conditions are set
  dt0 = 1; % the first index for array dt of time intervals
  dtf = nt - 1; % the last index for array dt of time intervals
  i0 = 1; % index for array t of time instances, beginning of interval dt(i) 
  i1 = 2; % index for array t of time instances, ending of interval dt(i)
otherwise
  error('d must be 1 or -1.')
end
% geting the dimentions of vectors and their initial values
if isempty(ic)
  ny = 0;
  y = zeros(nx,ny,nt); % declare the array of state variables
else
  if size(ic,1) ~= nx % check the dimention of initial conditions 
    error('IC must have the same number of rows as XMESH.')
  end
  ny = size(ic,2); % get number of state variables from IC second dimention
  y = zeros(nx,ny,nt); % declare the array of state variables
  y(:,:,i0) = ic; % write the initial conditions for all x
end

nz = size(icz,1); % get number of state variables from IC first dimention
z = zeros(nz,nt); % declare the array of state variables
if nz ~= 0
  z(:,i0) = icz; % write the initial conditions of z
  if d == 1, xb = zeros(nt,ny); end % kept to the next function call when d == -1
else xb_t = NaN;
end

if nargin > 8 % exist(u)
  if isempty(u), u = zeros(nx,0,nt); end
  if size(u,1) ~= nx, error('size(u,1) = nx'), end
else
  u = zeros(nx,0,nt);
end

if nargin > 9 % exist(v)
  if isempty(v), v = zeros(0,nt); end
else
  v = zeros(0,nt); 
end

if nargin > 7 %exist('fpq')
  [p0,q0] = feval(fpq,xmesh,t(i0),y(:,:,i0),z(:,i0),u(:,:,i0),v(:,i0),varargin{:});
  if isempty(p0) 
    np = 0;
    p0 = zeros(nx,0);
  else
    np = size(p0,2);
    if size(p0,1) ~= nx, error('size(p0,1) = nx'), end
  end
  if isempty(q0), nq = 0; 
  else nq = size(q0,1);
  end
end
p = zeros(nx,np,nt); % may be declare this variable as persistent ?????
q = zeros(nq,nt);  % may be declare this variable as persistent ?????
p(:,:,i0) = p0;
if nq ~= 0, q(:,i0) = q0;end

% arrays declaration before the loop for faster perfomance
xymesh = xmesh * ones(1,ny); % prepare rectangular mesh for state variables
y_int = zeros(nx,ny); % local array for integration (on rectangular mesh)
y_t = zeros(nx,ny);% local array of values on nonrectangular mesh
z_t = zeros(nz,1);% local array of z(t)
fz = z_t;% slops of characteristics
xymesh_t = zeros(nx,ny);% local array of nonrectangular mesh
% local array of boundary conditions
bc_t = zeros(ny);
if isnumeric(fbc)
  if isempty(fbc) % if there is no boundary conditions
    bc_t(:) = NaN; 
  else
    if length(fbc) ~= ny, error('length(fbc) = ny'), end
    bc_t = fbc;
  end
end
f = zeros(nx,ny);% right hand side values
c = ones(nx,ny);% slops of characteristics
% the main loop over time. nt-1 times
for i = dt0:d:dtf % index for array dt of time increments
  if d == 1, dti = dt(i); else dti = -dt(i); end % set sign to the time interval
  [f,c,fz] = feval(fpde,xmesh,t(i0),y(:,:,i0),z(:,i0),p(:,:,i0),q(:,i0),u(:,:,i0),v(:,i0),varargin{:}); % get functions f and c
  y_t = y(:,:,i0) + dti * f; % get approximations of y which generally are not on the grid nodes
  if nz ~= 0, z_t = z(:,i0) + dti * fz; end% get approximations of z
  if isempty(c), y_int = y_t;
  else
    if ~isnumeric(fbc) % get the boundary conditions for t(i1) but with q for t(i0) !
      bc_t = feval(fbc,t(i1),z_t,q(:,i0),v(:,i1),varargin{:});% and first order approx. of z
    end
    xymesh_t = xymesh + c * dti;% approximate x position of chracteristics on the next step t
%     if d == 1 && nz ~= 0, xb(i1,:) = 1-1./z_t; end % coordinates of discontinuety (persistent) ???????????????
    for j = 1:ny % quick interpolation columnwise on the xmesh
      if nz ~= 0, xb_t = xb(i1,j); end
      if xymesh(x0,j) == xymesh_t(x0,j) || isnan(bc_t(j)) || isinf(bc_t(j))
        y_int(:,j) = Interpolation(xymesh_t(:,j),y_t(:,j),xmesh,xb_t);
      else % interpolation with the use of boundary conditions
        if d == 1
          y_int(:,j) = Interpolation([xymesh(x0,j);xymesh_t(:,j)],[bc_t(j);y_t(:,j)],xmesh,xb_t);
        else
          y_int(:,j) = Interpolation([xymesh_t(:,j);xymesh(x0,j)],[y_t(:,j);bc_t(j)],xmesh,xb_t);
        end
      end
    end
  end
  % p and q
  if np ~= 0 || nq ~= 0 
    [p0,q0] = feval(fpq,xmesh,t(i1),y_int,z_t,u(:,:,i1),v(:,i1),varargin{:}); 
    if np ~= 0, p(:,:,i1) = p0;end
    if nq ~= 0, q(:,i1) = q0;end
  end
  % get functions f and c on the end of interval dt(i)
  [f,c,fz] = feval(fpde,xmesh,t(i1),y_int,z_t,p(:,:,i1),q(:,i1),u(:,:,i1),v(:,i1),varargin{:});
  if ~isempty(c)
    for j = 1:ny % quick interpolation columnwise on the xymesh_t
      f(:,j) = Interpolation(xmesh,f(:,j),xymesh_t(:,j),xb);
      if size(c,1) > 1, c(:,j) = Interpolation(xmesh,c(:,j),xymesh_t(:,j),xb);end % if c differs for x
    end
  end
  y_t = 0.5 * (y_t + y(:,:,i0) + dti * f); % take average between two approximations of y
  if nz ~= 0, z(:,i1) = 0.5 * (z_t + z(:,i0) + dti * fz); end% take average between two approximations of z
  if isempty(c), y(:,:,i1) = y_t;
  else
    if ~isnumeric(fbc)% get the boundary conditions for t(i1) and with q for t(i1) !
      bc_t = feval(fbc,t(i1),z(:,i1),q(:,i1),v(:,i1),varargin{:});
    end
    xymesh_t = 0.5 * (xymesh_t + xymesh + c * dti);% take average between two approximations of x positions
%     if d == 1 && nz ~= 0, xb(i1) = 1-1./z(:,i1);end % ???????????????
    for j = 1:ny % quick interpolation columnwise on the xmesh
      if nz ~= 0, xb_t = xb(i1,j); end
      if xymesh(x0,j) == xymesh_t(x0,j) || isnan(bc_t(j)) || isinf(bc_t(j))
        y(:,j,i1) = Interpolation(xymesh_t(:,j),y_t(:,j),xmesh,xb_t);
      else % interpolation with the use of boundary conditions
        if d == 1
          y(:,j,i1) = Interpolation([xymesh(x0,j);xymesh_t(:,j)],[bc_t(j);y_t(:,j)],xmesh,xb_t);
        else
          y(:,j,i1) = Interpolation([xymesh_t(:,j);xymesh(x0,j)],[y_t(:,j);bc_t(j)],xmesh,xb_t);
        end
      end
    end
  end
  % nonlocal variables p and q
  if np ~= 0 || nq ~= 0 
    [p0,q0] = feval(fpq,xmesh,t(i1),y(:,:,i1),z(:,i1),u(:,:,i1),v(:,i1),varargin{:}); 
    if np ~= 0, p(:,:,i1) = p0;end
    if nq ~= 0, q(:,i1) = q0;end
  end
  % get the boundary condition again with improved q 
  % it does not influence the order of convergance but increases accuracy
  % if boundary conditions depend on q
  if ~isnumeric(fbc) && ~isempty(c)
    bc_t = feval(fbc,t(i1),z(:,i1),q(:,i1),v(:,i1),varargin{:});
    for j = 1:ny % write the boundary conditions for t(it)
      if ~isnan(bc_t(j)) && ~isinf(bc_t(j)), y(x0,j,i1) = bc_t(j);end 
    end
  end
  %
  i0 = i0 + d;     % index for array t of time instances, beginning of interval dt(i) 
  i1 = i1 + d;     % index for array t of time instances, end of interval dt(i) 
end
end
function yi = Interpolation(x,y,xi,xb)
    % xp are the points of discontinuity (breaks)
    yi = pchip(x,y,xi);
    return
%     if any(isnan(x))
%         x(isnan(x)) = 0;
%     end
yi = zeros(size(xi));
method = 'linear'; % "!!!!!
if isnan(xb) || isempty(xb)
  yi = interp1(x,y,xi,method,'extrap');
else
  b = x <= xb;
  bi = xi <= xb;
  switch sum(b)
      case 0
      yi(bi) = y(1);
      case 1
      yi(bi) = y(b);
      otherwise        
      yi(bi) = interp1(x(b),y(b),xi(bi),method,'extrap');
  end
  switch sum(~b)
      case 0
      yi(~bi) = y(end);
      case 1
      yi(~bi) = y(~b);
      otherwise        
      yi(~bi) = interp1(x(~b),y(~b),xi(~bi),method,'extrap');
  end
end
end