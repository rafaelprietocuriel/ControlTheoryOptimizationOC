function display(mapPrim)

fixpointflag=isfixpoint(mapPrim);
structmapPrim=struct(mapPrim.doctrajectory);
try
    eigval=eig(mapPrim);
catch
    eigval=NaN;
end
if isequal(get(0,'FormatSpacing'),'compact')
    disp([inputname(1) ' =']);
    if ~isempty(mapPrim)
        disp('ocmatclass: mapprimitive')
        disp(['modelname: ' structmapPrim.modelname])
        if fixpointflag
            disp('Fixpoint:')
        else
            disp('Periodic Solution:')
        end
        if fixpointflag
            disp(structmapPrim.y(:,1));
        else
            disp(structmapPrim.y);
        end
        disp('Arcidentifier:')
        disp(structmapPrim.arcarg);
    else
        disp('empty mapprimitive');
    end
else
    disp(' ')
    disp([inputname(1) ' =']);
    disp(' ')
    if ~isempty(mapPrim)
        disp('ocmatclass: mapprimitive')
        disp(' ');
        disp(['modelname: ' structmapPrim.modelname])
        if fixpointflag
            disp(' ');
            disp('Fixpoint:')
        else
            disp(' ');
            disp('Period:')
            disp(mapPrim.period)
            disp(' ');
            disp('Periodic Solution:')
        end
        if fixpointflag
            disp(structmapPrim.y(:,1));
        else
            disp(structmapPrim.y);
        end
        disp(' ');
        disp('Eigenvalues:')
        if isempty(eigval)
            disp('  []')
        else
            disp(eigval);
        end
        disp(' ');
        disp('Arcidentifier:')
        disp(structmapPrim.arcarg);
    else
        disp('empty mapprimitive');
        disp(' ');
    end
end