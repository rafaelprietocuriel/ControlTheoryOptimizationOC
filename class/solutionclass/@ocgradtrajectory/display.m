function display(ocgTrj)

for ii=1:multiplicity(ocgTrj)
    if multiplicity(ocgTrj)>1
        disp(['Path ' num2str(ii)])
    end
    if isequal(get(0,'FormatSpacing'),'compact')
        disp([inputname(1) ' =']);
        if ~isempty(ocgTrj(ii))
            disp('ocmatclass: ocgradtrajectory')
            disp(struct(ocgTrj(ii)));
        else
            disp('empty ocgradtrajectory');
        end
    else
        disp(' ')
        disp([inputname(1) ' =']);
        disp(' ')
        if ~isempty(ocgTrj(ii))
            disp('ocmatclass: ocgradtrajectory')
            disp(' ');
            disp(struct(ocgTrj(ii)));
        else
            disp('empty ocgradtrajectory');
            disp(' ');
        end
    end
end