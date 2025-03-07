function [stateval meshidx]=state(ppdeObj,solObj,varargin)
%
% STATE returns the state values/variables.
%
% X=STATE(PPDEOBJ,SOLOBJ) SOLOBJ is a solution object
% (dynprimitive,octrajectory,ocasymptotic,occurve) a structure, including
% at least the fieldnames 'y' and 'arcposition'. If SOLOBJ is empty the
% output is the same as for STATE(PPDEOBJ). Otherwise the state values are
% returned. If SOLOBJ is an octrajectory consisting of multiple arcs X is a
% cell array of matrices, with the state values for each arc separately.

stateval=[];
meshidx=[];

if isempty(ppdeObj) || isempty(solObj)
    return
end
if ispdeprimitive(solObj) ||  ispdetrajectory(solObj) || ispdeasymptotic(solObj)
    femdat=femdata(solObj);
    arcarg=arcargument(solObj);
    indepvar=time(ppdeObj,solObj,1);
    depvar=dependentvar(solObj);
    arcpos=arcposition(solObj);
    arcn=arcnum(solObj);
end

stateval=zeros(ynum*npt,size(depvar,2));
for ii=1:size(depvar,2)
    for jj=1:ynum
        stateval(npt*(jj-1)+1:npt*jj,ii)=depvar(npt*(jj-1)+1:npt*jj,ii);
    end
end

