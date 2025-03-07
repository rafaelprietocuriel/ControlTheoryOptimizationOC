function display(pdePrim)

try
    eigval=eig(pdePrim);
catch
    eigval=NaN;
end
N=gridnum(pdePrim);
if isequal(get(0,'FormatSpacing'),'compact')
    disp([inputname(1) ' =']);
    if ~isempty(pdePrim)
        disp('ocmatclass: pdeprimitive')
        if isequilibrium(pdePrim)
            disp(['arcarg : ' int2str(pdePrim.pdetrajectory.arcarg)])
            disp('y : ')
            disp(reshape(pdePrim.pdetrajectory.y,N,[]).')
        else
            disp(pdePrim.pdetrajectory);
            disp('period')
            disp(pdePrim.period)
        end
        %disp(structppdePrim.arcarg);
    else
        disp('empty pdeprimitive');
    end
else
    disp(' ')
    disp([inputname(1) ' =']);
    disp(' ')
    if ~isempty(pdePrim)
        disp('ocmatclass: pdeprimitive')
        disp(' ');
        if isequilibrium(pdePrim)
            disp(['arcarg : ' int2str(pdePrim.pdetrajectory.arcarg)])
            disp(' ');
            disp('y : ')
            disp(reshape(pdePrim.pdetrajectory.y,N,[]).')
        else
            disp(pdePrim.pdetrajectory);
            disp(' ');
            disp('period')
            disp(pdePrim.period)
        end
        disp(' ');
        disp('Eigenvalues:')
        if isempty(eigval)
            disp('  []')
        else
            ns=length(find(real(eigval)<-1e-6));
            nc=length(find(abs(real(eigval))<1e-6));
            nu=length(eigval)-ns-nc;
            disp(['stable   : ' num2str(ns)]);
            disp(['unstable : ' num2str(nu)]);
            disp(['center   : ' num2str(nc)]);
        end
        disp(' ');
        %disp('Arcidentifier:')
        %disp(structppdePrim.arcarg);
    else
        disp('empty pdeprimitive');
        disp(' ');
    end
end