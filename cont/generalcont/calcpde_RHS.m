function res=calcpde_RHS(tmesh,y,freepar,modelpar,ode,bc)

global OCMATCONT

switch OCMATCONT.bvpmethod
    case {'bvp4c'}
        res=calcrespde_bvp4c(tmesh,y,freepar,modelpar,ode,bc);
    case {'bvp5c'}
        res=calcres_bvp5c(tmesh,y,freepar,modelpar,ode,bc);
    case {'bvp6c'}
        res=calcrespde_bvp6c(tmesh,y,freepar,modelpar,ode,bc);
end
