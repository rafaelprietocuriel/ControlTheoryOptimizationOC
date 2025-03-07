function ocTrjMP=setoptimaltimeindex(ocTrjMP,optindex)

if isempty(ocTrjMP)
    return
end

ocTrjMP.optimaltransition=optindex;