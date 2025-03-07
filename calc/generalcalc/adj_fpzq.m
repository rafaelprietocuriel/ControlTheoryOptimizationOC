function [eta,zeta] = adj_fpzq(x,t,xi,theta,uyp,vzq,varargin)
persistent u y z p v q nx ny np nv nq aayq expr aaxyq aaxyp aaxqp
% get control w
w = varargin{1};
% get options
opt = varargin{2};
% get dimentions
nx = length(x);
nu = opt.Nu;
ny = opt.Ny;
nz = opt.Nz;
np = opt.Np;
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
zeta = zeros(nq,1); % start constructing zeta
if isfield(opt,'funEtaZeta')
    [eta,zeta] = feval(opt.funEtaZeta,x,t,y,z,p,q,xi,theta,u,v,w,varargin{2:end});
else
    if isfield(opt,'funPhi_q') % size(funPhi_q) = [ny,nq]
        if ~isempty(opt.funPhi_q)
            aayq = zeros(ny,nq);  % auxiliary array
            if isnumeric(opt.funPhi_q), aayq(:,:) = opt.funPhi_q;
            else aayq(:,:) = feval(opt.funPhi_q,x,t,y,z,p,q,u,v,varargin{2:end});
            end
            zeta = zeta + (xi(1,:) * aayq)'; % vector sum along ny for all q
        end
    end
    expr = zeros(nx,nq); % expression for integration
    if isfield(opt,'funL_q') % size(funL_q) = [nx,nq]
        if ~isempty(opt.funL_q)
            if isnumeric(opt.funL_q), expr = expr + opt.funL_q;
            else expr = expr + feval(opt.funL_q,x,t,y,z,p,q,u,v,w,varargin{2:end});
            end
        end
    end
    if isfield(opt,'funF_q') % size(funF_q) = [nx,ny,nq]
        if ~isempty(opt.funF_q)
            aaxyq = zeros(nx,ny,nq); % auxiliary array
            if isnumeric(opt.funF_q), aaxyq(:,:,:) = opt.funF_q;
            else aaxyq(:,:,:) = feval(opt.funF_q,x,t,y,z,p,q,u,v,varargin{2:end});
            end
            for i=1:nq  % for all q takes vector sum along ny
                expr(:,i) = expr(:,i) + sum(xi .* aaxyq(:,:,i), 2);
            end
        end
    end
    zeta = zeta + trapzint2q(diff(x),expr);  %trapezoidal integration
    %
    eta = zeros(nx,np); % start constructing eta
    if isfield(opt,'funL_p')% size(funL_p) = [nx,np]
        if ~isempty(opt.funL_p)
            if isnumeric(opt.funL_p), eta = eta + opt.funL_p;
            else eta = eta + feval(opt.funL_p,x,t,y,z,p,q,u,v,w,varargin{2:end});
            end
        end
    end
    if isfield(opt,'funF_p')% size(funF_p) = [nx,ny,np]
        if ~isempty(opt.funF_p)
            aaxyp = zeros(nx,ny,np); % auxiliary array
            if isnumeric(opt.funF_p), aaxyp(:,:,:) = opt.funF_p;
            else aaxyp(:,:,:) = feval(opt.funF_p,x,t,y,z,p,q,u,v,varargin{2:end});
            end
            for i=1:np % for all p takes vector sum along ny
                eta(:,i) = eta(:,i) + sum(xi .* aaxyp(:,:,i), 2);
            end
        end
    end
    if isfield(opt,'funH_p')% size(funH_p) = [nx,nq,np]
        if ~isempty(opt.funH_p)
            aaxqp = zeros(nx,nq,np); % auxiliary array
            if isnumeric(opt.funH_p), aaxqp(:,:,:) = opt.funH_p;
            else aaxqp(:,:,:) = feval(opt.funH_p,x,t,y,z,q,u,v,varargin{2:end});
            end
            for i=1:np % for all p takes vector sum along nq
                eta(:,i) = eta(:,i) + aaxqp(:,:,i) * zeta;
            end
        end
    end
end