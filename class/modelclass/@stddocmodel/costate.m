function varargout=costate(docObj,solObj,varargin)
%
% COSTATE returns the costate values/variables.
% 
% L=COSTATE(OCOBJ) OCOBJ is a stdocmodel class. L is a cell array of
% strings consisting of the costate variable names. 
% 
% L=COSTATE(OCOBJ,SOLOBJ) SOLOBJ is a solution object
% (dynprimitive,octrajectory,ocasymptotic,occurve) a structure, including 
% at least the fieldnames 'y' and 'arcposition'. If SOLOBJ is empty the
% output is the same as for COSTATE(OCOBJ). Otherwise the costate values
% are returned. If SOLOBJ is an octrajectory consisting of multiple arcs L
% is a cell array of matrices, with the costate values for each arc
% separately. 
%
% L=COSTATE(OCOBJ,SOLOBJ,CONNECF) CONNECF=0 is the same as
% COSTATE(OCOBJ,SOLOBJ). If CONNECF is 1 multiple arcs are returned
% "connected", i.e., the costate values of all arcs are returned in one
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
    depvar=dependentvariable(docObj,arcarg);
    varargout{1}=depvar(costatecoord(docObj));
else
    for ii=1:arcn
        arcp=arcpos(1,ii):arcpos(2,ii)+1;
        ycoord=costatecoord(docObj);
        if connectflag
            varargout{1}(xcoord,arcp(1:end-1))=depvar(ycoord,arcp(2:end));
        else
            varargout{ii}=depvar(ycoord,arcp(2:end));
        end
    end
end