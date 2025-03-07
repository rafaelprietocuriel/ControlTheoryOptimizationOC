function disp(docTrj)

if isequal(get(0,'FormatSpacing'),'compact')
    if ~isempty(docTrj)
        disp('ocmatclass: doctrajectory')
        disp(struct(docTrj));
    else
        disp('empty doctrajectory');
    end
else
    if ~isempty(docTrj)
        disp('ocmatclass: doctrajectory')
        disp(' ');
        disp(struct(docTrj));
    else
        disp('empty doctrajectory');
        disp(' ');
    end
end