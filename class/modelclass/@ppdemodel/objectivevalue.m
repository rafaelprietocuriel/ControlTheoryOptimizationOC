function o=objectivevalue(ppdeObj,solObj,arcarg,varargin)
%
% 
o=[];
recalculate=[];
addstationary=[];
if isempty(ppdeObj)
    return
end
if nargin<=2
    arcarg=0;
end
if nargin==1
    solObj=[];
end

% check input arguments
recalculateidx=find(strcmp(varargin,'recalculate'));
if ~isempty(recalculateidx)
    recalculate=varargin{recalculateidx+1};
end
if isempty(recalculate)
    recalculate=0;
end

addstationaryidx=find(strcmp(varargin,'addstationary'));
if ~isempty(addstationaryidx)
    addstationary=varargin{addstationaryidx+1};
end
if isempty(addstationary)
    addstationary=1;
end

par=parametervalue(ppdeObj);
if isppdeprimitive(solObj)
    arcarg=0;
    t=[0 1];
    depvar=[solObj.y solObj.y];
    pt=points(solObj);
    eg=edges(solObj);
    tr=triangles(solObj);
    trianglearea=pdetrg(pt,tr);
    epdeSol=solObj;
elseif isppdetrajectory(solObj) || isppdeasymptotic(solObj)
    arcarg=0;
    t=time(ppdeObj,solObj,1);
    depvar=genstates(solObj);
    pt=points(solObj);
    eg=edges(solObj);
    tr=triangles(solObj);
    trianglearea=pdetrg(pt,tr);
    if isppdeasymptotic(solObj)
        epdeSol=ppdeprimitive(solObj);
    else
        epdeSol=ppdeprimitive([]);
    end
end
if isppdeasymptotic(solObj) && isfield(solObj.ppdetrajectory.discretizationinfo.coord,'objectivevalue') && ~isempty(solObj.ppdetrajectory.discretizationinfo.coord.objectivevalue) && ~recalculate
    o=solObj.y(solObj.ppdetrajectory.discretizationinfo.coord.objectivevalue,end);
elseif ~isppdeprimitive(solObj) && isppdetrajectory(solObj) && isfield(solObj.discretizationinfo.coord,'objectivevalue') && ~isempty(solObj.ppdetrajectory.discretizationinfo.coord.objectivevalue) && ~recalculate
    o=solObj.y(solObj.discretizationinfo.coord.objectivevalue,end);
else
    n=length(t);
    of=zeros(1,n);
    for ii=1:n
        depvarint=pdeintrp(eg,tr,depvar(:,ii));
        tmp=feval(ppdeObj,'ObjectiveFunction',t(ii),pt,depvarint,par,arcarg);
        of(ii)=sum(trianglearea.*tmp);
    end
    o=sum(diff(t).*(of(1:n-1)+of(2:n))/2);
end

if ~isempty(epdeSol) && strcmp(discountfactor(ppdeObj),'expdisc') && addstationary
    r=discountrate(ppdeObj);
    depvarint=pdeintrp(eg,tr,depvar(:,end));
    tmp=feval(ppdeObj,'ObjectiveFunction',0,pt,depvarint,par,arcarg);
    if r>0
        o=o+exp(-r*t(end))*sum(trianglearea.*tmp)/r;
    end
end