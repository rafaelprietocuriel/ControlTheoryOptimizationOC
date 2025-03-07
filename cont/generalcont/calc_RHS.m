function res=calc_RHS(tmesh,y,z,freepar,modelpar,ode,bc,icfun)

global OCMATCONT

switch OCMATCONT.bvpmethod
    case {'bvp4c'}
        res=calcres_bvp4c(tmesh,y,freepar,modelpar,ode,bc,icfun);
    case {'gbvp4c'}
        res=calcres_gbvp4c(tmesh,y,freepar,modelpar,ode,bc,icfun);
    case {'bvp5c'}
        res=calcres_bvp5c(tmesh,y,freepar,modelpar,ode,bc,icfun);
    case {'bvp6c'}
        res=calcres_bvp6c(tmesh,y,freepar,modelpar,ode,bc,icfun);
    case 'sbvpoc'
        res=calcres_sbvpoc(tmesh,y,z,freepar,modelpar,ode,bc,icfun);
end
