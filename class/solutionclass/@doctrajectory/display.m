function display(docTrj)

if isequal(get(0,'FormatSpacing'),'compact')
    disp([inputname(1) ' =']);
    if ~isempty(docTrj)
        disp('ocmatclass: doctrajectory')
        disp(struct(docTrj));
    else
        disp('empty doctrajectory');
    end
else
    disp(' ')
    disp([inputname(1) ' =']);
    disp(' ')
    if ~isempty(docTrj)
        disp('ocmatclass: doctrajectory')
        disp(' ');
        disp(struct(docTrj));
    else
        disp('empty doctrajectory');
        disp(' ');
    end
end