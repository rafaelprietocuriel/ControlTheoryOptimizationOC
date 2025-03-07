function ocTrj=uminus(ocTrj)

if isempty(ocTrj)
    return
end
ocTrj.y=-ocTrj.y;