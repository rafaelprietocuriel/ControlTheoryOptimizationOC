function varargout=time(dgObj,solObj,varargin)
%

connectflag=[];
if isempty(dgObj)
    return
end
if nargin>=3
    connectflag=varargin{1};
end
if isempty(connectflag)
    connectflag=0;
end
if isstruct(solObj)
    try
        try
            indepvar=solObj.independentvar;
        catch
            indepvar=solObj.x*solObj.arcinterval(end);
        end
        try
            depvar=solObj.dependentvar;
        catch
            depvar=solObj.y;
        end
        arcpos=solObj.arcposition;
        arcint=solObj.arcinterval;
        try
            inftimetransform=solObj.inftimetransformation;
        catch
            inftimetransform=[];
        end
        try
            x0=solObj.initialtime;
        catch
            x0=solObj.x0;
        end
        arcn=numel(solObj.arcarg);
    catch
        ocmaterror(['Field is missing.\n' lasterr])
    end
elseif isoctrajectory(solObj) || isocasymptotic(solObj)
    indepvar=independentvar(solObj);
    arcpos=arcposition(solObj);
    arcint=arcinterval(solObj);
    inftimetransform=inftimetransformation(solObj);
    x0=initialtime(solObj);
    arcn=arcnum(solObj);
end

diffarcinterval=diff(arcint);
t=x0;
for ii=1:arcn
    arcp=arcpos(1,ii):arcpos(2,ii);
    if isempty(inftimetransform) || ~inftimetransform
        t=t(end)+(indepvar(1,arcp)-ii+1)*diffarcinterval(ii);
    else
    end
    if connectflag
        varargout{1}(1,arcp)=t;
    else
        varargout{ii}=t;
    end
end