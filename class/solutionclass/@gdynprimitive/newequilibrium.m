function [gdynPrim,X,newarcid]=newequilibrium(gdynPrim,opt)
%
% returns the state-costate jacobian of an equilibrium.
X=[];
newarcid=[];
if isequilibrium(gdynPrim)
    ocObj=loadmodel(gdynPrim);
    if nargin==1
        opt=defaultocoptions;
    end
    [b,infoS]=testadmissibility(gdynPrim,ocObj,opt);
    arcid=arcargument(gdynPrim);

    if b==0
        newarcid=arcid;
        X=dependentvar(gdynPrim);
        X=X{1};
        ocmatmsg('Equilibrium is admissible.\n')
        return
    end
    n=statenum(ocObj);
    violationidx=find(infoS.violationmat);
    controlcnum=controlconstraintnum(ocObj);
    ccidx=constraintcombinationindex(ocObj,arcid);
    if violationidx<=controlcnum
        ccidx(violationidx)=1;
    elseif violationidx<=2*controlcnum
        ccidx(violationidx-controlcnum)=0;
    else
        return
    end
    newarcid=constraintcombinationindex2arcarg(ocObj,ccidx);
    coord=implicitcontrolcoordinate(ocObj,newarcid);
    X=dependentvar(gdynPrim);
    X=X{1}(1:2*n);
    ctrl=control(gdynPrim);
    ctrl=ctrl{1};
    X=[X;ctrl(coord)];
    gdynPrim=generateequilibrium(X,newarcid,modelparameter(gdynPrim),modelname(gdynPrim),class(gdynPrim));
else
    ocmatmsg('Not yet implemented for limit-cycles.\n')
end
