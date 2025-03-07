function arcpos=getarcposition(tmesh)
persistent N arcposition

if isempty(N) || length(tmesh)~=N
    N=length(tmesh);
    idx=find(diff(tmesh)==0);
    arcposition=[1 idx+1;idx N];
end

arcpos=arcposition;
