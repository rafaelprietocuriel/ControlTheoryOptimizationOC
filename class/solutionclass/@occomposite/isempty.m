function b=isempty(ocComp)
%
%
b=all(cellfun('isempty',ocComp.path));