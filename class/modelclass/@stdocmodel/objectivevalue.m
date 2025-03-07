function o=objectivevalue(ocObj,solObj,varargin)
%
%
o=[];
indepvar=0;
if isempty(ocObj) || isempty(solObj)
    return
end
if nargin<=2
    arcarg=0;
end
if nargin==1
    solObj=[];
end
par=parametervalue(ocObj);
if isdynprimitive(solObj) 
    arcarg=arcargument(solObj);
    indepvar=time(ocObj,solObj,1);
    depvar=dependentvar(solObj);
    arcpos=arcposition(solObj);
    arcn=arcnum(solObj);
elseif  isocasymptotic(solObj) || isoctrajectory(solObj) 
    arcarg=arcargument(solObj);
    indepvar=time(ocObj,solObj,1);
    depvar=dependentvar(solObj);
    arcpos=arcposition(solObj);
    arcn=arcnum(solObj);
elseif isstruct(solObj) 
    arcarg=solObj.arcarg;
    indepvar=solObj.x*solObj.arcinterval(end);
    depvar=solObj.y;
    if isfield(solObj,'arcposition')
        arcpos=solObj.arcposition;
    else
        arcpos=[1;length(indepvar)];
    end
    arcn=length(arcarg);
end
o=0;

if isfield(solObj.solverinfo,'objectivevaluecoord')
    o=depvar(solObj.solverinfo.objectivevaluecoord,end);
elseif isautonomous(ocObj)
    if isocasymptotic(solObj)
        r=discountrate(ocObj);
        for ii=1:arcn
            arcp=arcpos(1,ii):arcpos(2,ii);
            o(arcp)=feval(ocObj,'Hamiltonian',indepvar(arcp),depvar(:,arcp),par,arcarg(ii))/r;
        end
    elseif isoctrajectory(solObj)
        o=feval(ocObj,'Objectivevalue',indepvar,depvar,par,arcarg);
    end
    if isoctrajectory(solObj) && ~isocasymptotic(solObj)
        try
            S=feval(ocObj,'DiscountedSalvagevalue',indepvar(end),depvar(:,end),par,arcarg(end));
            o=o+S;
        end
    end
else
    for ii=1:arcn
        arcp=arcpos(1,ii):arcpos(2,ii);
        of=feval(ocObj,'ObjectiveFunction',indepvar(arcp),depvar(:,arcp),par,arcarg(ii));
        %o=o+sum(diff(indepvar).*(of(1:end-1)+of(2:end))/2);
        o(arcp)=[o(end) o(end)+cumsum(diff(indepvar(arcp)).*(of(1:end-1)+of(2:end))/2)];
    end
    if isoctrajectory(solObj) && ~isocasymptotic(solObj)
        try
            S=feval(ocObj,'DiscountedSalvagevalue',indepvar(end),depvar(:,end),par,arcarg(end));
            o=o+S;
        end
    end
end