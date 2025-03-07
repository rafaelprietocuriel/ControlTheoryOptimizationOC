function b=isidentical(ocTrj1,ocTrj2)
%
%
b=0;
if isempty(ocTrj1) && isempty(ocTrj2)
    b=1;
   return
end

if ~isempty(ocTrj1) && ~isempty(ocTrj2)
    arcarg1=arcargument(ocTrj1);
    arcarg2=arcargument(ocTrj2);
    if length(arcarg1)~=length(arcarg2) || ~all(arcarg1==arcarg2)
        return
    end
    if ~all(arcinterval(ocTrj1)==arcinterval(ocTrj2))
        return
    end
    if ~all(all(arcposition(ocTrj1)==arcposition(ocTrj2)))
        return
    end
    x1=independentvar(ocTrj1);
    x2=independentvar(ocTrj2);
    if length(x1)~=length(x2) || ~all(x1==x2)
        return
    end
    y1=dependentvar(ocTrj1);
    y2=dependentvar(ocTrj2);
    if length(y1)~=length(y2) || ~all(all(y1==y2))
        return
    end
    b=1;
else
    return
end