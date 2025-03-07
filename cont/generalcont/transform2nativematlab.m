function [sol,freepar,n,N]=transform2nativematlab(tmesh,coeff,modelpar)
global OCMATCONT OCBVP
N=numel(tmesh);
sol=[];
freepar=[];
n=[];
switch OCMATCONT.bvpmethod
    case {'bvp4c','bvp6c'}
        n=OCBVP.numode;
        nN=n*N;
        y=reshape(coeff(1:nN),n,N);
        freepar=coeff(nN+(1:OCMATCONT.HE.numparameter));
        if strcmp(OCMATCONT.bvpmethod,'bvp4c')
            yp=calcdata_bvp4c(tmesh,y,freepar,modelpar,OCMATCONT.ode);
        elseif strcmp(OCMATCONT.bvpmethod,'bvp6c')
            [yp,ypmid]=calcdata_bvp6c(tmesh,y,freepar,modelpar,OCMATCONT.ode);
            sol.ypmid=ypmid;
        end
        sol.x=tmesh;
        sol.y=y;
        sol.yp=yp;
        sol.solver=OCMATCONT.bvpmethod;
        if isfield(OCMATCONT.HE,'numparametermc')
            sol.parameters=freepar(1:OCMATCONT.HE.numparametermc);
        end
    case {'gbvp4c'}
        y=zeros(OCBVP.maxnumode,length(tmesh));
        y(OCMATCONT.HE.DDATA.meshvalcoord)=coeff(OCMATCONT.HE.ycoord);
        freepar=coeff(OCMATCONT.HE.parametercoord); % parameter values
        %y=reshape(coeff(1:nN),n,N);
        %freepar=coeff(nN+(1:OCMATCONT.HE.numparameter));
        yp=calcdata_gbvp4c(tmesh,y,freepar,modelpar,OCMATCONT.ode);
        sol.x=tmesh;
        sol.y=y;
        sol.yp=yp;
        sol.solver=OCMATCONT.bvpmethod;
        if isfield(OCMATCONT.HE,'numparametermc')
            sol.parameters=freepar(1:OCMATCONT.HE.numparametermc);
        end
    case 'bvp5c'
        n=OCBVP.numode;
        neqn=OCBVP.neqn;
        y=reshape(coeff,neqn,N);
        paramcoeff=n+(1:OCMATCONT.HE.numparameter);
        freepar=coeff(paramcoeff);
        [yp ymid]=calcdata_bvp5c(tmesh,y,freepar,modelpar,OCMATCONT.ode);
        y(paramcoeff,:)=[];
        yp(paramcoeff,:)=[];
        ymid(paramcoeff,:)=[];
        sol.x=tmesh(1:OCBVP.nstages:end);
        sol.y=y(:,1:OCBVP.nstages:end);
        sol.solverinfo.yp=yp;
        sol.solverinfo.ymid=ymid;
        sol.solver=OCMATCONT.bvpmethod;

    case 'sbvpoc'
        sol.x=tmesh(OCMATCONT.MESHDATA.tmeshidx);
        sol.y=coeff2points(tmesh,coeff,'grid');
        N=OCMATCONT.MESHDATA.meshNumber;
        sol.parameters=coeff(OCMATCONT.MESHDATA.freeparameterindex);
    otherwise
        return
end

