function b=isidentical(ocgTrj1,ocgTrj2)
%
%
b=0;
if isempty(ocgTrj1) && isempty(ocgTrj2)
    b=1;
   return
end

if ~isempty(ocgTrj1) && ~isempty(ocgTrj2)
    arcarg1=arcargument(ocgTrj1);
    arcarg2=arcargument(ocgTrj2);
    if length(arcarg1)~=length(arcarg2) || ~all(arcarg1==arcarg2)
        return
    end
    if ~all(arcinterval(ocgTrj1)==arcinterval(ocgTrj2))
        return
    end
    if ~all(all(arcposition(ocgTrj1)==arcposition(ocgTrj2)))
        return
    end
    x1=independentvar(ocgTrj1);
    x2=independentvar(ocgTrj2);
    if length(x1)~=length(x2) || ~all(x1==x2)
        return
    end
    y1=dependentvar(ocgTrj1);
    y2=dependentvar(ocgTrj2);
    if length(y1)~=length(y2) 
        return
    end
    for ii=1:length(y1)
        if ~all(all(y1{ii}==y2{ii}))
            return
        end
    end
    b=1;
else
    return
end