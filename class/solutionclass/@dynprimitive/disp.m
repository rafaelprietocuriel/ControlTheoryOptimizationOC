function disp(dynPrim)

structdynPrim=struct(dynPrim.octrajectory);
try
    eigval=eig(dynPrim);
catch
    eigval=[];
end
if isequal(get(0,'FormatSpacing'),'compact')
    if ~isempty(dynPrim)
        disp('ocmatclass: dynprimitive')
        disp(['modelname: ' structdynPrim.modelname])
        disp(structdynPrim.y);
    else
        disp('empty dynprimitive');
    end
else
    if ~isempty(dynPrim)
        disp('ocmatclass: dynprimitive')
        disp(' ');
        disp(['modelname: ' structdynPrim.modelname])
        disp(' ');
        disp(structdynPrim.y);
        disp(' ');
        disp('Eigenvalues:')
        disp(eigval);
    else
        disp('empty dynprimitive');
        disp(' ');
    end
end