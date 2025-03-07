function disp(ocgTrj)

for ii=1:multiplicity(ocgTrj)
    if multiplicity(ocgTrj)>1
        disp(['Path ' num2str(ii)])
    end
    if isequal(get(0,'FormatSpacing'),'compact')
        if ~isempty(ocgTrj)
            disp('ocmatclass: ocgradtrajectory')
            disp(struct(ocgTrj));
        else
            disp('empty octrajectory');
        end
    else
        if ~isempty(ocgTrj)
            disp('ocmatclass: ocgradtrajectory')
            disp(' ');
            disp(struct(ocgTrj));
        else
            disp('empty ocgradtrajectory');
            disp(' ');
        end
    end
end