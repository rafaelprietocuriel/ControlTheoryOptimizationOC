function varargout=state(docObj,solObj,varargin)
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
if isempty(docObj)
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
if ismapprimitive(solObj) || isdoctrajectory(solObj) || isdocasymptotic(solObj)
    depvar=[initialstate(solObj) dependentvar(solObj)];
    arcpos=arcposition(solObj);
    arcn=arcnum(solObj);
elseif isoccurve(solObj)
    depvar=dependentvar(solObj);
    arcpos=[0;size(depvar,2)-1];
    arcn=1;
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
    % remove initial point
    depvar=dependentvariable(docObj,arcarg);
    varargout{1}=depvar(statecoord(docObj),2:end);
else
    for ii=1:arcn
        arcp=arcpos(1,ii):arcpos(2,ii)+1;
        ycoord=statecoord(docObj);
        if connectflag
            varargout{1}(ycoord,arcp(1:end-1))=depvar(ycoord,arcp(2:end));
        else
            varargout{ii}=depvar(ycoord,arcp(2:end));
        end
    end
end