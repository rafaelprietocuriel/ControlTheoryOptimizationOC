function [sol_coarse,sol_fine] = errorestimate4dae(sol_coarse)
global OCMATCONT

MESHDATA=OCMATCONT.MESHDATA;

[tcolmesh1_2,coeff1_2]=coeff2coeff1_2(sol_coarse.data.xcol,sol_coarse.data.coeff);

OCMATCONT.MESHDATA=OCMATCONT.MESHDATA1_2;
try
    % since coeff1_2 is computed on the fine mesh, the global variable
    % OCMATCONT has to be adapted, i.e. MESHDATA has to be replace by
    % MESHDATA1_2. After the computation the change has to be undone
    [tcolmesh1_2,coeff1_2]=newtoncorrection(tcolmesh1_2,coeff1_2,[]);
catch
    OCMATCONT.MESHDATA=MESHDATA;
    rethrow(lasterror)
end

OCMATCONT.MESHDATA=MESHDATA;
if isempty(coeff1_2)
    sol_coarse=[]; 
    sol_fine=[];
    return
end

%y1=coeff2points(sol_coarse.data.xcol,sol_coarse.data.coeff,'collocationgrid2collocationgrid1_2'); %coarse solution on fine grid + collpoints
y1=coeff2points(sol_coarse.data.xcol,sol_coarse.data.coeff,'general',tcolmesh1_2); %coarse solution on fine grid + collpoints

sol_fine=sol2daestruct(1,tcolmesh1_2,coeff1_2);
%y2_2=coeff2points(sol_fine.data.xcol,coeff1_2,'collocationgrid1_22collocationgrid'); %fine solution on coarse grid + collpoints
y2_2=coeff2points(sol_fine.data.xcol,coeff1_2,'general1_2',sol_coarse.data.xcol); %fine solution on coarse grid + collpoints

%fine solution on coarse grid + collpoints
sol_coarse.data.errest=(2^OCMATCONT.CollocationNumber/(1-2^OCMATCONT.CollocationNumber))*(y2_2-sol_coarse.data.ycol);

%coarse solution on fine grid + collpoints
sol_fine.data.errest=(1/(1-2^OCMATCONT.CollocationNumber))*(sol_fine.data.ycol-y1);

end
