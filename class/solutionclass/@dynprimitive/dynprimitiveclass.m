function classname=dynprimitiveclass(dynPrim)

classname='';
if isempty(dynPrim)
    return
end

if ~dynPrim.period && size(dynPrim.octrajectory.y,2)==1
    classname='equilibrium';
elseif ~dynPrim.period && size(dynPrim.octrajectory.y,2)==2
    classname='fixpoint';
else
end