function out=dynprimitive(ocgAsym)

if isempty(ocgAsym)
    out=dynprimitve();
    return
end

out=ocgAsym.dynprimitive;