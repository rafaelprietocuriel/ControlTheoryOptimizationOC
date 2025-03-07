function b=isperiodic(dynPrim)
%
%
b=~isempty(dynPrim.period)&&(dynPrim.period~=0);