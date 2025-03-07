function out=subsref(mmObj,S)

switch S(1).type
    case '()'
        out=mmObj.Model{S(1).subs{1}};
    case '.'
        out=mmObj.(S(1).subs);
end
