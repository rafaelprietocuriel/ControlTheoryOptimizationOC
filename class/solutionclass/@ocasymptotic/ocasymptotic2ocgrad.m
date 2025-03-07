function ocgAsym=ocasymptotic2ocgrad(ocAsym)

if isempty(ocAsym)
    ocgAsym=ocasymptotic2ocgrad([]);
    return
end

if isempty(modelname(ocAsym))
    ocgAsym=ocasymptotic2ocgrad([]);
    return
end

ocgEP=dynprimitive2ocgrad(limitset(ocAsym));
ocgAsym=ocgradasymptotic(octrajectory2ocgrad(octrajectory(ocAsym)),ocgEP);