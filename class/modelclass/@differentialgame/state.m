function varargout=state(dgObj,solObj,varargin)
%
% STATE returns the state values/variables.
% 
% X=STATE(OCOBJ) OCOBJ is a stdocmodel class. X is a cell array of strings
% consisting of the state variable names. 
% 
% X=STATE(OCOBJ,SOLOBJ) SOLOBJ is a solution object
% (dynprimitive,octrajectory,ocasymptotic,occurve) a structure, including 
% at least the fieldnames 'y' and 'arcposition'. If SOLOBJ is empty the
% output is the same as for STATE(OCOBJ). Otherwise the state values are
% returned. If SOLOBJ is an octrajectory consisting of multiple arcs X is a
% cell array of matrices, with the state values for each arc separately.
%
% X=STATE(OCOBJ,SOLOBJ,CONNECF) CONNECF=0 is the same as
% STATE(OCOBJ,SOLOBJ). If CONNECF is 1 multiple arcs are returned
% "connected", i.e., the state values of all arcs are returned in one
% matrix.
connectflag=[];
if isempty(dgObj)
    return
end
if nargin==1
    solObj=[];
end
if nargin>=3
    connectflag=varargin{1};
end
if isempty(connectflag)
    connectflag=0;
end
par=parametervalue(dgObj);
if isstruct(solObj)
    try
        depvar=solObj.y;
        arcpos=solObj.arcposition;
    catch
        ocmaterror('If second input argument is a structure the fieldnames ''y'' and ''arcposition'' have to exist!')
    end
elseif isdynprimitive(solObj) || isoctrajectory(solObj) || isocasymptotic(solObj)
    arcarg=arcargument(solObj);
    indepvar=time(dgObj,solObj,1);
    depvar=dependentvar(solObj);
    arcpos=arcposition(solObj);
    arcn=arcnum(solObj);
elseif isoccurve(solObj)
    arcn=1;
    depvar=dependentvar(solObj);
    arcpos=arcposition(solObj);
    if isempty(arcpos)
        arcpos=[1 size(depvar,2)]';
    end
end
if isempty(solObj)
    % return symbolic expressions
    % if repflag=1 control values/Lagrange multipliers are replaced by the
    % corresponding symbolic expressions for the optimal values
    arcarg=[];
    if nargin>=3
        arcarg=varargin{1};
    end
    if isempty(arcarg)
        arcarg=0;
    end
    depvar=dependentvariable(dgObj,arcarg);
    varargout{1}=depvar(statecoord(dgObj));
else
    for ii=1:arcn
        arcp=arcpos(1,ii):arcpos(2,ii);
        xcoord=statecoordinate(dgObj);
        if 1%~stateconstraintflag
            x=depvar(xcoord,arcp);
        else
            x=feval(dgObj,'State',indepvar(arcp),depvar(:,arcp),par,arcarg(ii));
        end
        if connectflag
            varargout{1}(xcoord,arcp)=x;
        else
            varargout{ii}=x;
        end
    end
end