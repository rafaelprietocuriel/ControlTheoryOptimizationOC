function sol=addaecomponent(dynPrim,ocObj,sol,solvername,initnummesh,domaindata)

if isempty(dynPrim)
    return
end
if isequilibrium(dynPrim)
    domaindata=domaindata(arcargument(dynPrim)+1);
    if isempty(domaindata.implicitcontrolindex)
        return
    end
    switch solvername
        case 'sbvpoc'
            ctrl=control(ocObj,dynPrim);
            sol.y(domaindata.implicitcontrolcoord,1:initnummesh)=repmat(ctrl(domaindata.implicitcontrolindex),1,initnummesh);
        otherwise
    end
end