function disp(gdynPrim)

structgdynPrim=struct(gdynPrim.ocgtrajectory);
try
    eigval=eig(gdynPrim);
catch
    eigval=[];
end
if isequal(get(0,'FormatSpacing'),'compact')
    if ~isempty(gdynPrim)
        disp('ocmatclass: gdynprimitive')
        disp(['modelname: ' structgdynPrim.octrajectory.modelname])
        disp(structgdynPrim.octrajectory.y{1});
        disp('OdeNum:')
        disp(structgdynPrim.odenum);
    else
        disp('empty gdynprimitive');
    end
else
    if ~isempty(gdynPrim)
        disp('ocmatclass: gdynprimitive')
        disp(' ');
        disp(['modelname: ' structgdynPrim.octrajectory.modelname])
        disp(' ');
        disp(structgdynPrim.octrajectory.y{1});
        disp('OdeNum:')
        disp(' ');
        disp(structgdynPrim.odenum);
        disp(' ');
        disp('Eigenvalues:')
        disp(eigval);
    else
        disp('empty dynprimitive');
        disp(' ');
    end
end