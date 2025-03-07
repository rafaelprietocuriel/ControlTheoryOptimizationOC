function [odefinal,bcfinal,jacfinal,bcjacfinal,Joptions,dBCoptions] = ...
    bvpcontsetting(odefun,bcfun,odejac,bcjac,tmesh,coeff,tangent)
%BVPFUNCTIONS  Helper function for processing functional arguments for BVP solvers.
%
%   See also BVP4C, BVP5C, BVPSET.

%   Copyright 2007 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $  $Date: 2007/05/23 18:54:21 $

global OCMATCONT OCBVP
bvpmethod=OCMATCONT.bvpmethod;
OCBVP.nregions
odefinal = ode;
bcfinal = bc;
jacfinal = odejac;
bcjacfinal = bcjac;

% incorporate known parameters (for BVP4C)
if ~isempty(extras) && strcmp(bvpmethod,'bvp4c')
    if nparam == 0
        if nregions == 1
            odefinal = @(x,y) ode(x,y,extras{:});
            if isa(odejac,'function_handle')
                jacfinal = @(x,y) odejac(x,y,extras{:});
            end
        else
            odefinal = @(x,y,region) ode(x,y,region,extras{:});
            if isa(odejac,'function_handle')
                jacfinal = @(x,y,region) odejac(x,y,region,extras{:});
            end
        end
        bcfinal = @(ya,yb) bc(ya,yb,extras{:});
        if isa(bcjac,'function_handle')
            bcjacfinal = @(ya,yb) bcjac(ya,yb,extras{:});
        end
    else
        if nregions == 1
            odefinal = @(x,y,p) ode(x,y,p,extras{:});
            if isa(odejac,'function_handle')
                jacfinal = @(x,y,p) odejac(x,y,p,extras{:});
            end
        else
            odefinal = @(x,y,region,p) ode(x,y,region,p,extras{:});
            if isa(odejac,'function_handle')
                jacfinal = @(x,y,region,p) odejac(x,y,region,p,extras{:});
            end
        end
        bcfinal = @(ya,yb,p) bc(ya,yb,p,extras{:});
        if isa(bcjac,'function_handle')
            bcjacfinal = @(ya,yb,p) bcjac(ya,yb,p,extras{:});
        end
    end
end

% incorporate unknown parameters (for BVP5C)
if (nparam > 0) && strcmp(bvpmethod,'bvp5c')
    if nregions == 1
        odefinal = @odeParameters;
        if ~isempty(odejac)
            jacfinal = @jacParameters;
        end
    else
        msg = sprintf('Multipoint BVP not yet supported with BVP5C. Use BVP4C instead.');
        error('MATLAB:bvpfunctions:MultipointBVPNotSupported','%s',msg);
    end
    bcfinal  = @bcParameters;
    if ~isempty(bcjac)
        bcjacfinal = @bcjacParameters;
    end
end

% data for numerical Jacobians
switch bvpmethod
    case 'bvp4c'
        total_eqns = neqn;
    case 'bvp5c'
        total_eqns = neqn + nparam;  % add equations for unknown parameters
end
Joptions = [];
dBCoptions = [];
threshval = 1e-6;
if isempty(odejac)
    Joptions.diffvar = 2;  % dF(x,y)/dy
    vectorized = strcmp(bvpget(options,'Vectorized','off'),'on');
    if vectorized  % xy-vectorized
        Joptions.vectvars = [1,2];
    else
        Joptions.vectvars = [];
    end
    Joptions.thresh = threshval(ones(total_eqns,1));
end
if isempty(bcjac)
    dBCoptions.vectvars = [];
    dBCoptions.thresh = threshval(ones(total_eqns,1));
    dBCoptions.fac_dya = [];
    dBCoptions.fac_dyb = [];
end

% ---------------------------------------------------------
% Nested functions
% ---------------------------------------------------------

    function f = odeParameters(x,y)
        p = y(neqn+1:neqn+nparam,1);  % extract p from y(:,1)
        % add trivial equations for unknown parameters
        f = [ ode(x,y(1:neqn,:),p);
            zeros(nparam,numel(x))];
    end  % odeParameters

% ---------------------------------------------------------

    function res = bcParameters(ya,yb)
        p = ya(neqn+1:neqn+nparam);   % extract p from ya
        res = bc(ya(1:neqn),yb(1:neqn),p);
    end  % bcParameters

% ---------------------------------------------------------

    function J = jacParameters(x,y)
        if iscell(odejac)  % constant Jacobians
            dfdy = odejac{1};
            dfdp = odejac{2};
        else
            p = y(neqn+1:neqn+nparam,1);   % extract p from y(:,1)
            [dfdy,dfdp] = odejac(x,y(1:neqn,:),p);
        end
        % add trivial equations for unknown parameters
        J = [dfdy, dfdp;
            zeros(nparam,neqn+nparam)];
    end  % jacParameters

% ---------------------------------------------------------

    function [dya,dyb] = bcjacParameters(ya,yb)
        if iscell(bcjac)  % constant Jacobians
            dbcdya = bcjac{1};
            dbcdyb = bcjac{2};
            dbcdp  = bcjac{3};
        else
            p = ya(neqn+1:neqn+nparam);   % extract p from ya
            [dbcdya,dbcdyb,dbcdp] = bcjac(ya(1:neqn,:),yb(1:neqn,:),p);
        end
        dya = [dbcdya, dbcdp];
        dyb = [dbcdyb, dbcdp];
    end  % bcjacParameters

% ---------------------------------------------------------

end  % bvpfunctions
