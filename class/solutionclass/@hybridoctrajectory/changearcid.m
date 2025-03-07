function hocTrj=changearcid(hocTrj,changearcidx,newarcid)

if isempty(hocTrj) || isempty(changearcidx)
    return
end

if length(changearcidx)~=length(newarcid)
    return
end
hocTrj.arcarg(changearcidx)=newarcid;