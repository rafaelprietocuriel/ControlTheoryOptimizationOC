function Jac=calcpde_RHSJac(tmesh,y,freepar,modelpar,ode,bc,odejac,bcjac)

global OCMATCONT

switch OCMATCONT.bvpmethod
    case {'bvp4c'}
        Jac=calcresjacpde_bvp4c(tmesh,y,freepar,modelpar,ode,bc,odejac,bcjac);
    case {'bvp5c'}
        Jac=calcresjac_bvp5c(tmesh,y,freepar,modelpar,ode,bc,[],odejac,bcjac,[]);
    case {'bvp6c'}
        Jac=calcresjacpde_bvp6c(tmesh,y,freepar,modelpar,ode,bc,odejac,bcjac);
end
