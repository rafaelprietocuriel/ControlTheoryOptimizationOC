function cellstring=string2cell(string,matrixtype)
% each cell for matrixtype 'matrix' consists of a row of the matrix
% each cell for matrixtype 'vector' consists of an entry of the vector

cellstring=[];
string=removematrixstring(string);
if ~strcmp(string([1 end]),'[]')
    cellstring{1}=string;
    return
end
switch matrixtype
    case 'vector'
        cellstring=replacediff(string(2:end-1));
        cellstring=regexp(cellstring,',','split');
    case 'matrix'
        cellstring=replacediff(string(2:end-1));
        cellstring=regexp(cellstring,'],(\ )*[','split');
    case 'charmatrix'
        cellstring=replacediff(string(2:end-1));
        cellstring=regexp(cellstring,'],.[','split');
end

function str=replacediff(str)

while 1
    S=regexp(str, 'diff\(','once');
    if isempty(S)
        break
    end
    open  = 0;
    start = S+5;
    for kk = start : numel(str)
        switch str(kk)
            case '('
                open = open+1;
            case ')'
                if open
                    open = open-1;
                else
                    match=str(start:kk-1);
                    arg=regexp(match,',','split','once');
                    S1=regexp(arg{1}, '\(','once');
                    if ~isempty(S1)
                        open1  = false;
                        start1 = S1+1;
                        for mm = start1 : numel(arg{1})
                            switch arg{1}(mm)
                                case '('
                                    open1 = true;
                                case ')'
                                    if open1
                                        open1 = false;
                                    else
                                        funcname=arg{1}(1:start1-2);
                                        S2=regexp(arg{2},'\$','once');
                                        if S2
                                            tmp=arg{2};
                                            while S2
                                                S3=regexp(tmp,'\)','once');
                                                arg3=regexp(tmp(S2+3:S3-1),',','split','once');
                                                reptmp='';
                                                for jj=1:str2num(arg3{2})
                                                    reptmp=[reptmp arg3{1} ','];
                                                end
                                                arg{2}=strrep(arg{2},tmp(S2-1:S3),reptmp(1:end-1));
                                                S2=regexp(arg{2},'\$','once');
                                            end
                                        end
                                        derarg=regexp(arg{2},',','split');
                                        if isempty(derarg)
                                            derarg=arg(2);
                                        end
                                        replacestr=funcname;
                                        for ii=1:length(derarg)
                                            replacestr=[ replacestr '_D' derarg{ii}];
                                        end
                                        replacestr=['D' num2str(length(derarg)) replacestr  arg{1}(start1-1:mm)];
                                        break;
                                    end
                            end
                        end
                    end
                    break;
                end
        end
    end
    str=strrep(str,['diff(' match ')'],replacestr);
    
end