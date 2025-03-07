function disp(ocgMTrj)

if isequal(get(0,'FormatSpacing'),'compact')
    if ~isempty(ocgMTrj)
        disp('ocmatclass: ocgradmultipath')
        disp(['degree   : ' int2str(degree(ocgMTrj))])
        disp(struct(ocgMTrj));
    else
        disp('empty ocgradmultipath');
    end
else
    if ~isempty(ocgMTrj)
        disp('ocmatclass: ocgradmultipath')
        disp(' ');
        disp(['degree   : ' int2str(degree(ocgMTrj))])
        disp(' ');
        disp(struct(ocgMTrj));
    else
        disp('empty ocgradmultipath');
        disp(' ');
    end
end