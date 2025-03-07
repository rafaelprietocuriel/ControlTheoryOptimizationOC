function [y,z] = goded(d,fode,ic,icz,xmesh,t,u,v,varargin)
%   [Y,Z] = GNLQLODE1(d,FODE,IC,XMESH,TSPAN,PFUN) solves
%   initial value problems for small systems of local quasi-linear ODEs 
%   in one space variable x and time t  
%   simultaniously with local ODEs in time t.
%
%   Dy(x,t)/Dt = f(x,t,y(x,t),z(t))
%
%   where vector f(x,t,y(x,t),z(t)) 
%   can depend on nonlocal functions with respect to space x
%
%   Dz(t)/Dt = fz(x,t,y(x,t),z(t)) 
%
%   function fz can contain any integral trasformations along
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
%   [FY,FZ] = FODE(XMESH,T,Y,Z,U,V) is a function that evaluates the quantities
%   defining the differential equation. The input arguments are
%   XMESH is the vector of values x
%   T is the scalat time
%   Y is the matrix of state variables. size(Y) = [nx,ny]
%   Z is the vector. size(Z) = [nz,1]
%   FODE returns matrices
%   FY is the right hand side of the nonlocal ODE. size(FY) = [nx,ny]
%   FZ is the right hand side of the local ODE. size(FZ) = [nz,1]
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
%   GNLQLODE1 returns values of the solution on a mesh provided as the input
%   array XMESH. The entries of XMESH must satisfy 
%       a = XMESH(1) < XMESH(2) < ... < XMESH(NX) = b 
%   for some NX >= 3.  The entries of TSPAN must satisfy 
%       t0 = TSPAN(1) < TSPAN(2) < ... < TSPAN(NT) = tf 
%   for some NT >= 3. 
%  
%   GNLQLDE1 returns matrices
%   Y is matrix of state variables y, size(Y) = [nx,ny,nt]
%   Z is matrix of variable z, size(Z) = [nz,nt]
%
%   GNLQLODE1 uses Heun's method (also called the modified Euler's method)
%   to solve ODE. 
%
%   Vladimir M. Veliov
%   Anton O. Belyakov
%   Revision date: 2010/10/25

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
if any(dt <= 0)
  error('The entries of TSPAN must be strictly increasing.')
end

%xmesh = xmesh(:); % make shure that xmesh is a one-dimentional array
nx = length(xmesh);
if nx < 3
  error('XMESH must have at least 3 entries.')
end
h = diff(xmesh); % argument increments
if any(h <= 0)
  error('The entries of XMESH must be strictly increasing.')
end
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
z(:,i0) = icz; % write the initial conditions of z

if nargin > 6 % exist(u)
  if isempty(u), u = zeros(nx,0,nt); end
  if size(u,1) ~= nx, error('size(u,1) = nx'), end
else u = zeros(nx,0,nt);
end

if nargin > 7 % exist(v)
  if isempty(v), v = zeros(0,nt); end
else
  v = zeros(0,nt); 
end

% arrays declaration before the loop for faster perfomance
y_t = zeros(nx,ny);% local array of values y(t) on mesh
z_t = zeros(nz);% local array of z(t)
fy = y_t;% right hand side values
fz = z_t;% slops of characteristics
% the main loop over time. nt-1 times
for i = dt0:d:dtf % index for array dt of time increments
  if d == 1, dti = dt(i); else dti = -dt(i); end % set sign to the time interval
  [fy,fz] = feval(fode,xmesh,t(i0),y(:,:,i0),z(:,i0),u(:,:,i0),v(:,i0),varargin{:}); % get functions f and c
  y_t = y(:,:,i0) + dti * fy; % get approximations of y
  z_t = z(:,i0) + dti * fz; % get approximations of z
  % get functions fy and fz on the end of interval dt(i)
  [fy,fz] = feval(fode,xmesh,t(i1),y_t,z_t,u(:,:,i1),v(:,i1),varargin{:});
  y(:,:,i1) = 0.5 * (y_t + y(:,:,i0) + dti * fy); % take average between two approximations of y
  z(:,i1) = 0.5 * (z_t + z(:,i0) + dti * fz); % take average between two approximations of z
  %
  i0 = i0 + d;     % index for array t of time instances, beginning of interval dt(i) 
  i1 = i1 + d;     % index for array t of time instances, ending of interval dt(i) 
end