function b=isimplicit(dgObj,arcid)
%
% ISIMPLICIT test if implicit control exist.
%
% ISIMPLICIT(dgObj) test if for any arcid any control is implicit.
%
% ISIMPLICIT(dgObj,ARCID) if ARCID is empty it works like
% ISIMPLICIT(dgObj). Otherwise it tests if any control for the specified
% ARCID is implicit.
%

b=false;
if isempty(dgObj)
    return
end

if nargin==1
    arcid=[];
end

if isempty(arcid)
    arcid=0;
    while 1
        if isimplicit(dgObj,arcid)
            b=true;
            break
        end
        arcid=arcid+1;
        if arcid==arcnum(dgObj)
            break
        end
    end
else
    if any(implicitcontrolcoordinate(dgObj,arcid))
        b=true;
    end
end