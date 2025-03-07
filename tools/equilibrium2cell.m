function varargout=equilibrium2cell(ocEP)

if iscell(ocEP)
    return
end
if isa(ocEP,'gdynprimitive') || isa(ocEP,'dynprimitive')
    ocEP0{1}=ocEP;
    assignin('caller',inputname(1),ocEP0);
    if nargout==1
        varargout{1}=true;
    end
end