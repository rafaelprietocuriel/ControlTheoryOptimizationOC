function out=discountedsalvagevalue(ocObj,solObj,varargin)
%
% SALVAGEVALUE

if isempty(ocObj)
    out=[];
    return
end
if nargin==1
    solObj=[];
end
par=parametervalue(ocObj);
if isstruct(solObj)
    try
        arcarg=solObj.arcarg;
        indepvar=solObj.independentvar;
        depvar=solObj.dependentvar;
    catch
        ocmaterror('If the second input argument is a structure the fields ''dependentvar'', ''independentvar'' and ''arcarg'' have to exist!')
    end
elseif isdaeoctrajectory(solObj)
    arcarg=[];
    indepvar=independentvar(solObj);
    depvar=dependentvar(solObj);
elseif isdynprimitive(solObj) || isoctrajectory(solObj) || isocasymptotic(solObj)    
    arcarg=arcargument(solObj);
    indepvar=time(ocObj,solObj,1);
    depvar=dependentvar(solObj);
    arcpos=arcposition(solObj);
    arcn=arcnum(solObj);
elseif  isoccurve(solObj)
    arcarg=arcargument(solObj);
    indepvar=zeros(size(arcarg));
    depvar=dependentvar(solObj);
    arcpos=arcposition(solObj);
    if isempty(arcpos)
        arcarg=unique(arcarg);
        arcpos=[1 size(depvar,2)]';
        arcn=1;
    else
        arcpos=arcposition(solObj);
        arcn=arcnum(solObj);
    end

end

if isempty(solObj)
    try
        varargout{1}=feval(ocObj,'DiscountedSalvagevalue',[],[],par);
    catch
        varargout{1}=[];
        return
    end
else
    % return value of the canonical system evaluated at 'depvar'
    out=feval(ocObj,'DiscountedSalvagevalue',indepvar(end),depvar(:,end),par,arcarg);
end
