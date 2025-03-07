function varargout=state(ocObj,solObj,varargin)
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
if isempty(ocObj)
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
xcoord=statecoord(ocObj);
if ischar(connectflag)
    depvar=jumpdependentvar(solObj);
    jumparg=jumpargument(solObj);
    jumparg([1 end])=-1;
    switch lower(connectflag(1))
        case 'l' % left side limit
            depvar=depvar(:,1:2:2*length(jumparg)-1);
        case 'r' % right side limit
            depvar=depvar(:,2:2:2*length(jumparg));
    end
    depvar(:,jumparg==0)=[];
    varargout{1}=depvar(xcoord,:);
else
    if isstruct(solObj)
        try
            depvar=solObj.y;
            arcpos=solObj.arcposition;
        catch
            ocmaterror('If second input argument is a structure the fieldnames ''y'' and ''arcposition'' have to exist!')
        end
    elseif ishybridoctrajectory(solObj)
        arcarg=arcargument(solObj);
        indepvar=time(ocObj,solObj,1);
        depvar=dependentvar(solObj);
        arcpos=arcposition(solObj);
        %arcpos(1,1)=1;
        %arcpos(2,end)=arcpos(2,end)+1;
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
        depvar=dependentvariable(ocObj,arcarg);
        varargout{1}=depvar(statecoord(ocObj));
    else
        for ii=1:arcn
            arcp=arcpos(1,ii):arcpos(2,ii);
            x=depvar(xcoord,arcp);
            if connectflag
                varargout{1}(xcoord,arcp)=x;
            else
                varargout{ii}=x;
            end
        end
    end
end