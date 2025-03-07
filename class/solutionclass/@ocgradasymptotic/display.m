function display(ocgAsym)

for ii=1:multiplicity(ocgAsym)
    if multiplicity(ocgAsym)>1
        disp(['Path ' num2str(ii)])
    end
    if isequal(get(0,'FormatSpacing'),'compact')
        disp([inputname(1) ' =']);
        if ~isempty(ocgAsym(ii))
            disp('ocmatclass: ocgradasymptotic')
            disp(ocgAsym(ii));
        else
            disp('empty ocgradasymptotic');
        end
    else
        disp(' ')
        disp([inputname(1) ' =']);
        disp(' ')
        if ~isempty(ocgAsym(ii))
            disp('ocmatclass: ocgradasymptotic')
            disp(' ');
            disp(ocgAsym(ii).limitset);
            disp(ocgAsym(ii).ocgradtrajectory);
        else
            disp('empty ocgradasymptotic');
            disp(' ');
        end
    end
end