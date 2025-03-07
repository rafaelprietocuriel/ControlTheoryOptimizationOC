function Jac=calc_monodromy(tmesh,y,z,freepar,modelpar,ode,odejac)

global OCMATCONT

switch OCMATCONT.bvpmethod
    case {'bvp4c'}
        [F Fmid]=calcdata_bvp4c(tmesh,y,freepar,modelpar,ode);
        Jac=calcmonodromy_bvp4c(tmesh,y,freepar,modelpar,ode,odejac,F,Fmid);
    case {'bvp5c'}
        Jac=calcmonodromy_bvp5c(tmesh,y,freepar,modelpar,ode,odejac);
    case {'bvp6c'}
        [F Fmid Xmid Ymid]=calcdata_bvp6c(tmesh,y,freepar,modelpar,ode);
        Jac=calcmonodromy_bvp6c(tmesh,y,freepar,modelpar,ode,odejac,F,Fmid,Xmid,Ymid);
    case 'sbvpoc'
        boundarystate(tmesh,y,z);
        Jac=calcmonodromy_sbvpoc(tmesh,y,z,freepar,modelpar,ode,odejac);
end
