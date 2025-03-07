function [u,v,w,Val,du,dv,dw,z,theta] = fun_gradmap_gode(u,v,w,gmopt,varargin)
persistent aazv
%persistent z theta %ht
Val = NaN;
% get options
opt = varargin{1};
% get initial conditions for primal system
% get initial conditions for primal system z
if isfield(opt,'funZ0')
    if ~isempty(opt.funZ0)
        if isnumeric(opt.funZ0)
            icz = opt.funZ0;
        else
            icz = feval(opt.funZ0,w,varargin{:});
        end
    else error('vector funZ0 must be nonempty')
    end
else error('function or vector funZ0 must be given')
end
% solve primal system
if isfield(opt,'funF')
    z = gode(1,opt.funF,icz,opt.tspan,v,w,varargin{:});
else
    z = gode(1,opt.funS,icz,opt.tspan,v,varargin{:});
end
if isfield(opt,'MixConstr') % mixed constraints to control depending on state variables
    [u,v,w] = feval(opt.MixConstr,opt.tspan,z,u,v,w,varargin{:});
end
% calculate value function
if isfield(gmopt,'CalcVal')
    if gmopt.CalcVal
        %if isfield(opt,'dt'), ht = opt.dt;
        %else ht = diff(opt.tspan); end
        if isfield(opt,'funJ')
            Val = feval(opt.funJ,opt.tspan,z,v,w,varargin{:});
        else
            Val = 0;
            if isfield(opt,'funLT')
                if ~isempty(opt.funLT)
                    Val = Val + feval(opt.funLT,z(:,end),varargin{:});
                end
            end
            if isfield(opt,'funLt')
                if ~isempty(opt.funLt)
                    Val = Val + trapzint1q(opt.dt,feval(opt.funLt,opt.tspan,z,v,w,varargin{:}));
                end
            end
        end
    end
end
% galculate gradient of the hamiltonian
if isfield(gmopt,'CalcGrad')
    if gmopt.CalcGrad == 0
        %du = Val;
        du = [];
        dv = [];
        dw = [];
        return
    end
    % by default calculate gradient % else return
end
% get initial conditions for adjoint system
if isfield(opt,'funl_z')
    if ~isempty(opt.funl_z)
        if isnumeric(opt.funl_z)
            icz = opt.funl_z;
        else
            icz = feval(opt.funl_z,z(:,end),v(:,end),w,varargin{:});
        end
    else error('vector funl_z must be nonempty')
    end
else error('function or vector funl_z must be given')
end
% solve adjoint system regarding u and concatenated v,z along dimention 1.
% control w is passed as an additional parameter
theta = gode(-1,@adj_fode,icz,opt.tspan,cat(1,v,z),w,varargin{:});
% calculates the gradient of Hamiltonians
if isfield(opt,'funGradHam')
    [du,dv,dw] = feval(opt.funGradHam,opt.tspan,z,theta,v,w,varargin{:});
else
    % constract the gradiant of Hamiltonians by derivatives of the functions
    nv = opt.Nv;
    nz = opt.Nz;
    nt = opt.Nt;
    dv = zeros(nv,nt);
    if isfield(opt,'funS_v')
        if ~isempty(opt.funS_v)
            if isnumeric(opt.funS_v)
                if ndims(opt.funS_v) == 2, dv(:,:) = (theta' * opt.funS_v)';
                else % if ndims(opt.funS_v) == 3 i.e. size(opt.funS_v) = [nz,nv,nt]
                    for i = 1:nv
                        dv(i,:) = sum(theta .* reshape(opt.funS_v(:,i,:),[nz,nt]),1);
                    end
                end
            else
                aazv = zeros(nz,nv); % auxiliary array of Jacobian matrix
                for j = 1:nt
                    aazv(:,:) = feval(opt.funS_v,opt.tspan(j),z(:,j),v(:,j),varargin{:});
                    dv(:,j) = dv(:,j) + aazv' * theta(:,j); % Jacobian matrix is transposed !
                end
            end
        end
    end
    if isfield(opt,'funL_v')
        if isnumeric(opt.funL_v)
            if size(opt.funL_v,2) == nt, dv = dv + opt.funL_v;
            else % if size(opt.funL_v,2) == 1 i.e. size(opt.funL_v) = [nv,1]
                for i = 1:nv
                    dv(i,:) = dv(i,:) +  opt.funL_v(i);
                end
            end
        else
            for j = 1:nt
                dv(:,j) = dv(:,j) + feval(opt.funL_v,opt.tspan(j),z(:,j),v(:,j),w,varargin{:});
            end
        end
    end
    du = [];
    dw = [];
end