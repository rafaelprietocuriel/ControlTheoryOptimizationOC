function Jac=calc_RHSJac(tmesh,y,z,freepar,modelpar,ode,bc,icfun,odejac,bcjac,icfunjac)

global OCMATCONT

switch OCMATCONT.bvpmethod
    case {'bvp4c'}
        Jac=calcresjac_bvp4c(tmesh,y,freepar,modelpar,ode,bc,icfun,odejac,bcjac,icfunjac);
    case {'gbvp4c'}
        Jac=calcresjac_gbvp4c(tmesh,y,freepar,modelpar,ode,bc,icfun,odejac,bcjac,icfunjac);
    case {'bvp5c'}
        Jac=calcresjac_bvp5c(tmesh,y,freepar,modelpar,ode,bc,icfun,odejac,bcjac,icfunjac);
    case {'bvp6c'}
        Jac=calcresjac_bvp6c(tmesh,y,freepar,modelpar,ode,bc,icfun,odejac,bcjac,icfunjac);
    case 'sbvpoc'
        boundarystate(tmesh,y,z);
        Jac=calcresjac_sbvpoc(tmesh,y,z,freepar,modelpar,ode,bc,icfun,odejac,bcjac,icfunjac);
end
