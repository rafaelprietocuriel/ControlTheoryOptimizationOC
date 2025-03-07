function b=isimplicit(ocObj,arcid)
%
% ISIMPLICIT test if implicit control exist.
%
% ISIMPLICIT(OCOBJ) test if for any arcid any control is implicit.
%
% ISIMPLICIT(OCOBJ,ARCID) if ARCID is empty it works like
% ISIMPLICIT(OCOBJ). Otherwise it tests if any control for the specified
% ARCID is implicit.
%

b=false;
if isempty(ocObj)
    return
end

if nargin==1
    arcid=[];
end

if isempty(arcid)
    arcid=0;
    while 1
        if isimplicit(ocObj,arcid)
            b=true;
            break
        end
        arcid=arcid+1;
        if arcid==arcnum(ocObj)
            break
        end
    end
else
    if any(implicitcontrolcoordinate(ocObj,arcid))
        b=true;
    end
end