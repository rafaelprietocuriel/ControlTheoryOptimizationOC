function ocTrj=transform2multiarc(ocTrj,varargin)

solverinfoStruct=solverinfo(ocTrj);

% multiarccalc 0 potential multi point BVP is transformed to teo pohint
% BVP, arcposition denotes rows for each arc
if isfield(solverinfoStruct,'multiarccalc')
    multiarccalc=solverinfoStruct.multiarccalc;
else
    multiarccalc=1;
end

if multiarccalc
    return
else
    depvar=dependentvar(ocTrj);
    t=independentvar(ocTrj);
    numarc=arcnum(ocTrj);
    arcpos=arcposition(ocTrj);
    numarccoord=1:numarc;
    ocTrj.x=(t(ones(numarc,1)))+numarccoord(:,ones(1,numel(t)));
    ocTrj.x=reshape(ocTrj.x.',1,[]);
    sizedepvar=size(depvar,2);
    ocTrj.arcposition=[(0:numarc-1)*sizedepvar+1; ...
        (1:numarc)*sizedepvar];
    counter=0;
    ocTrj.y=zeros(max(diff(arcpos))+1,numarc*sizedepvar);
    for ii=1:numarc
        counter_start=counter+1;
        counter=counter+sizedepvar;
        ocTrj.y(:,counter_start:counter)=depvar(arcpos(1,ii):arcpos(2,ii),:);
    end
    ocTrj.solverinfo=[];
    ocTrj.solverinfo.coeff=[];
    ocTrj.solverinfo.tangent=[];
    ocTrj.solverinfo.tmesh=[];
    ocTrj.solverinfo.tmesh.pathtype=solverinfoStruct.pathtype;
    ocTrj.solverinfo.inftimetransformation=solverinfoStruct.inftimetransformation;
    try
        ocTrj.solverinfo.parameters=solverinfoStruct.parameters;
    end
    try
        ocTrj.solverinfo.continuationparameter=solverinfoStruct.continuationparameter;
    end
    ocTrj.solverinfo.multiarccalc=1;
    ocTrj.solver=[];
end