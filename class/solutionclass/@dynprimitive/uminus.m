function dynPrim=uminus(dynPrim)

if isempty(dynPrim)
    return
end
dynPrim.y=-dynPrim.y;