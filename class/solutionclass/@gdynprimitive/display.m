function display(gdynPrim)

equilibflag=isequilibrium(gdynPrim);
structgdynPrim=struct(gdynPrim.ocgtrajectory);
try
    eigval=eig(gdynPrim);
catch
    eigval=NaN;
end
if isequal(get(0,'FormatSpacing'),'compact')
    disp([inputname(1) ' =']);
    if ~isempty(gdynPrim)
        disp('ocmatclass: gdynprimitive')
        disp(['modelname: ' structgdynPrim.octrajectory.modelname])
        if equilibflag
            disp('Equilibrium:')
            disp(structgdynPrim.octrajectory.y{1});
            disp('Arcidentifier:')
            disp(structgdynPrim.octrajectory.arcarg);
            disp('OdeNum:')
            disp(structgdynPrim.odenum);
        else
            disp('Periodic Solution:')
        end
    else
        disp('empty dynprimitive');
    end
else
    disp(' ')
    disp([inputname(1) ' =']);
    disp(' ')
    if ~isempty(gdynPrim)
        disp('ocmatclass: gdynprimitive')
        disp(' ');
        disp(['modelname: ' structgdynPrim.octrajectory.modelname])
        if equilibflag
            disp(' ');
            disp('Equilibrium:')
            disp(structgdynPrim.octrajectory.y{1});
            disp(' ');
            disp('Eigenvalues:')
            if isempty(eigval)
                disp('  []')
            else
                disp(eigval);
            end
            disp(' ');
            disp('Arcidentifier:')
            disp(structgdynPrim.octrajectory.arcarg);
            disp('OdeNum:')
            disp(structgdynPrim.odenum);
        else
            disp(' ');
            disp('Period:')
            disp(gdynPrim.period)
            disp(' ');
            disp('Periodic Solution:')
        end
    else
        disp('empty dynprimitive');
        disp(' ');
    end
end