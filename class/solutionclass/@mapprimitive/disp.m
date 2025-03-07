function disp(mapPrim)

structmapPrim=struct(mapPrim.doctrajectory);
try
    eigval=eig(mapPrim);
catch
    eigval=[];
end
if isequal(get(0,'FormatSpacing'),'compact')
    if ~isempty(mapPrim)
        disp('ocmatclass: mapprimitive')
        disp(['modelname: ' structmapPrim.modelname])
        disp(structmapPrim.y);
    else
        disp('empty mapprimitive');
    end
else
    if ~isempty(mapPrim)
        disp('ocmatclass: mapprimitive')
        disp(' ');
        disp(['modelname: ' structmapPrim.modelname])
        disp(' ');
        disp(structmapPrim.y);
        disp(' ');
        disp('Eigenvalues:')
        disp(eigval);
    else
        disp('empty mapprimitive');
        disp(' ');
    end
end