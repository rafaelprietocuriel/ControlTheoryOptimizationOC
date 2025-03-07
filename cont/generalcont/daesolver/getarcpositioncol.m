function arcpos=getarcpositioncol(tcolmesh)
persistent N arcposition

if isempty(N) || length(tcolmesh)~=N
    N=length(tcolmesh);
    idx=find(diff(tcolmesh)==0);
    arcposition=[1 idx+1;idx N];
end

arcpos=arcposition;
