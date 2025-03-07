function hocTrj=removejump(hocTrj,jumpidx,jumpid)

if isempty(jumpidx)
    hocTrj.jumparg(end)=jumpid;
else
    hocTrj.jumparg(jumpidx)=jumpid;
end