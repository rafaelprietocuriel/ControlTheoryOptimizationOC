function disp(statEx)

for ii=1:multiplicity(statEx)
    if multiplicity(statEx)>1
        disp(['Solution ' num2str(ii)])
    end
    if isequal(get(0,'FormatSpacing'),'compact')
        if ~isempty(statEx)
            disp('ocmatclass: staticextremal')
            disp(struct(statEx));
        else
            disp('empty staticextremal');
        end
    else
        if ~isempty(statEx)
            disp('ocmatclass: staticextremal')
            disp(' ');
            disp(struct(statEx));
        else
            disp('empty staticextremal');
            disp(' ');
        end
    end
end