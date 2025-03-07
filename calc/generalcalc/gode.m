function z = gode(d,fode,icz,t,v,varargin)
%   Z = GODE(d,FODE,ICZ,TSPAN) solves
%   initial value problems for small systems of local ODEs
%   in time t
%
%   Dz(t)/Dt = fz(x,t,y(x,t),z(t),q(t))
%
%   functions fz and fq can contain any integral trasformations along
%   the space x
%
%   Parameter d must be -1 or 1. It sets the direction of the problem
%   (positive or negative) simultaniously both in time t and argument x.
%
%   All multi-dimentional arrays have the following structure rule:
%   The first dimention is state variable index if any.
%   The last dimention is time index if any.
%
%   FZ = FODE(T,Z,V) is a function that evaluates the quantities
%   defining the differential equation. The input arguments are
%   T is the scalat time
%   Z is the vector. size(Z) = [nz,1]
%   FODE returns matrix:
%   FZ is the right hand side of the local ODE. size(FZ) = [nz,1]
%
%   ICZ is the vector containing the initial conditions. size(ICZ) = [nz,1]
%
%   TSPAN must be strictly increasing array of time inctances at which the
%   solution is evaluetad
%
%   GODE returns matrices
%   Z is matrix of variable z, size(Z) = [nz,nt]
%
%   GODE uses Heun's method (also called the modified Euler's method)
%   to solve ODE.
%
%   Vladimir M. Veliov
%   Anton O. Belyakov
%   Revision date: 2010/02/12
persistent nz fz z_t nt dt dt0 dtf i0 i1
% check input parameters
if nargin < 4
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
% define initial and final indices of TSPAN and dt
switch d % direction of the problem forward (d=1) or backward (d=-1)
    case -1
        dt0 = nt - 1; % the first index for array dt of time intervals
        dtf = 1; % the last index for array dt of time intervals
        i0 = nt; % index for array t of time instances, beginning of interval dt(i)
        i1 = nt - 1; % index for array t of time instances, ending of interval dt(i)
    case 1
        dt0 = 1; % the first index for array dt of time intervals
        dtf = nt - 1; % the last index for array dt of time intervals
        i0 = 1; % index for array t of time instances, beginning of interval dt(i)
        i1 = 2; % index for array t of time instances, ending of interval dt(i)
    otherwise
        error('d must be 1 or -1.')
end
% geting the dimentions of vectors and their initial values
nz = size(icz,1); % get number of state variables from IC first dimention
z = zeros(nz,nt); % declare the array of state variables
z(:,i0) = icz; % write the initial conditions of z

if nargin > 4 % exist(v)
    if isempty(v)
        v = zeros(0,nt); 
    end
else
    v = zeros(0,nt);
end

% arrays declaration before the loop for faster perfomance
z_t = zeros(nz);% local array of z(t)
fz = z_t;% slops of characteristics
% the main loop over time. nt-1 times
for i = dt0:d:dtf % index for array dt of time increments
    if d == 1
        dti = dt(i); 
    else
        dti = -dt(i); 
    end % set sign to the time interval
    fz = feval(fode,t(i0),z(:,i0),v(:,i0),varargin{:}); % get functions f and c
    z_t = z(:,i0) + dti * fz; % get approximations of z
    % get functions fy and fz on the end of interval dt(i)
    fz = feval(fode,t(i1),z_t,v(:,i1),varargin{:});
    z(:,i1) = 0.5 * (z_t + z(:,i0) + dti * fz); % take average between two approximations of z
    % index for array t of time instances
    i0 = i0 + d; % beginning of interval dt(i)
    i1 = i1 + d; % ending of interval dt(i)
end