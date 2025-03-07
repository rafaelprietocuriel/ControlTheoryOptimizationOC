function display(dynPrim)

equilibflag=isequilibrium(dynPrim);
structdynPrim=struct(dynPrim.octrajectory);
try
    eigval=eig(dynPrim);
catch
    eigval=NaN;
end
if isequal(get(0,'FormatSpacing'),'compact')
    disp([inputname(1) ' =']);
    if ~isempty(dynPrim)
        disp('ocmatclass: dynprimitive')
        disp(['modelname: ' structdynPrim.modelname])
        if equilibflag
            disp('Equilibrium:')
        else
            disp('Periodic Solution:')
        end
        disp(structdynPrim.y);
        disp('Arcidentifier:')
        disp(structdynPrim.arcarg);
    else
        disp('empty dynprimitive');
    end
else
    disp(' ')
    disp([inputname(1) ' =']);
    disp(' ')
    if ~isempty(dynPrim)
        disp('ocmatclass: dynprimitive')
        disp(' ');
        disp(['modelname: ' structdynPrim.modelname])
        if equilibflag
            disp(' ');
            disp('Equilibrium:')
        else
            disp(' ');
            disp('Period:')
            disp(dynPrim.period)
            disp(' ');
            disp('Periodic Solution:')
        end
        disp(structdynPrim.y);
        disp(' ');
        disp('Eigenvalues:')
        if isempty(eigval)
            disp('  []')
        else
            disp(eigval);
        end
        disp(' ');
        disp('Arcidentifier:')
        disp(structdynPrim.arcarg);
    else
        disp('empty dynprimitive');
        disp(' ');
    end
end