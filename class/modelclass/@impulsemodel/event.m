function out=event(ocObj,solObj,varargin)
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
if isempty(ocObj)
    return
end
if nargin==1
    solObj=[];
end
xcoord=statecoord(ocObj);
depvar=jumpdependentvar(solObj);
jumparg=jumpargument(solObj);
jumparg([1 end])=-1;
depvarl=depvar(xcoord,1:2:2*length(jumparg)-1);
depvarr=depvar(xcoord,2:2:2*length(jumparg));
depvarl(:,jumparg==0)=[];
depvarr(:,jumparg==0)=[];
out=depvarr-depvarl;
