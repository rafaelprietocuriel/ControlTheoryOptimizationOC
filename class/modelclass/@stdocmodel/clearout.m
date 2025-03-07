function clearout(ocObj,force)


if nargin==1
    force=0;
end
clearout(loadmodeldata(ocObj),force);