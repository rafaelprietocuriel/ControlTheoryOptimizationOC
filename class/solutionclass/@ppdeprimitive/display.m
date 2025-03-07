function display(ppdePrim)

try
    eigval=eig(ppdePrim);
catch
    eigval=NaN;
end

if isequal(get(0,'FormatSpacing'),'compact')
    disp([inputname(1) ' =']);
    if ~isempty(ppdePrim)
        disp('ocmatclass: ppdeprimitive')
        disp(ppdePrim.ppdetrajectory);
        disp('period')
        disp(ppdePrim.period)
        %disp(structppdePrim.arcarg);
    else
        disp('empty ppdeprimitive');
    end
else
    disp(' ')
    disp([inputname(1) ' =']);
    disp(' ')
    if ~isempty(ppdePrim)
        disp('ocmatclass: ppdeprimitive')
        disp(' ');
        disp(ppdePrim.ppdetrajectory);
        disp(' ');
        disp('period')
        disp(ppdePrim.period)
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
        disp('empty ppdeprimitive');
        disp(' ');
    end
end