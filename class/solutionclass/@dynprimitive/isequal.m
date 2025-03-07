function b=isequal(dynPrim1,dynPrim2)
%
%
b=false;
if isempty(dynPrim1)
    if isempty(dynPrim2)
        b=true;
    end
    return
elseif isempty(dynPrim2)
    return
end
if isequilibrium(dynPrim1)
    if isequilibrium(dynPrim2)
        b=length(dynPrim1.octrajectory.y)==length(dynPrim2.octrajectory.y) && all(dynPrim1.octrajectory.y-dynPrim2.octrajectory.y==0) && dynPrim1.octrajectory.arcarg-dynPrim2.octrajectory.arcarg==0; % ...
            %&& strcmp(dynPrim1.octrajectory.modelname,dynPrim2.octrajectory.modelname) && all(dynPrim1.octrajectory.modelparameter-dynPrim2.octrajectory.modelparameter==0 | isnan(dynPrim1.octrajectory.modelparameter-dynPrim2.octrajectory.modelparameter));
    else 
        return
    end
end
if isequilibrium(dynPrim2)
    return
end
