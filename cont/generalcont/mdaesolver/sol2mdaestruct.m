function sol=sol2daestruct(flag,tcolmesh,coeff)
global OCMATCONT
switch flag
    case 0
        tmesh=tcolmesh(OCMATCONT.MESHDATA.tmeshidx);

        sol.x=tmesh;
        sol.y=coeff2points(tcolmesh,coeff,'grid');
        sol.parameters=coeff(OCMATCONT.MESHDATA.freeparameterindex);
        sol.data.xcol=tcolmesh;
        sol.data.ycol=coeff2points(tcolmesh,coeff,'collocationgrid');
        sol.data.coeff=coeff;
    case 1
        tmesh=tcolmesh(OCMATCONT.MESHDATA1_2.tmeshidx);
        sol.x=tmesh;
        sol.y=coeff2points(tcolmesh,coeff,'grid1_2');
        sol.parameters=coeff(OCMATCONT.MESHDATA1_2.freeparameterindex);
        sol.data.xcol=tcolmesh;
        sol.data.ycol=coeff2points(tcolmesh,coeff,'general1_2',tcolmesh);
        sol.data.coeff=coeff;
end
