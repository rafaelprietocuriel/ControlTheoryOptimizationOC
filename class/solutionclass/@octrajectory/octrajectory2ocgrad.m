function ocgTrj=octrajectory2ocgrad(ocTrj)

if isempty(ocTrj)
    ocgTrj=ocgradtrajectoy([]);
    return
end

if isempty(modelname(ocTrj))
    ocgTrj=ocgradtrajectoy([]);
    return
end
ocObj=stdocmodel(modelname(ocTrj));
par=modelparameter(ocTrj);
ocObj=changeparametervalue(ocObj,par);

ocgTrj=struct(ocTrj);
ocgTrj.t=time(ocObj,ocTrj,1);
[ocgTrj.t,idx]=unique(ocgTrj.t);
ocgTrj.y=state(ocObj,ocTrj,1);
ocgTrj.y=ocgTrj.y(:,idx);
ocgTrj.cst_y=costate(ocObj,ocTrj,1);
ocgTrj.cst_y=ocgTrj.cst_y(:,idx);
ocgTrj.v=control(ocObj,ocTrj,1);
ocgTrj.v=ocgTrj.v(:,idx);
ocgTrj.objectivevalue=[];%objectivevalue(ocObj,ocTrj,1);
ocgTrj.timehorizon=timehorizon(ocTrj);
ocgTrj.t0=ocgTrj.x0;
ocgTrj=rmfield(ocgTrj,{'x0','x'});
ocgTrj=ocgradtrajectory(ocgTrj);