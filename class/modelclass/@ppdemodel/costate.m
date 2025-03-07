function costateval=costate(ppdeObj,solObj,varargin)
%
% STATE returns the state values/variables.
%
% X=STATE(PPDEOBJ,SOLOBJ) SOLOBJ is a solution object
% (dynprimitive,octrajectory,ocasymptotic,occurve) a structure, including
% at least the fieldnames 'y' and 'arcposition'. If SOLOBJ is empty the
% output is the same as for STATE(PPDEOBJ). Otherwise the state values are
% returned. If SOLOBJ is an octrajectory consisting of multiple arcs X is a
% cell array of matrices, with the state values for each arc separately.

costateval=[];

if isempty(ppdeObj) || isempty(solObj)
    return
end

if isppdeprimitive(solObj) ||  isppdetrajectory(solObj) || isppdeasymptotic(solObj)
    depvar=dependentvar(solObj);
    pt=points(solObj);
    npt=numpoints(solObj);
    t=time(ppdeObj,solObj);
    nt=length(t);
end
coeffidx=meshindex(ppdeObj,solObj,'costate');
costateval=zeros(statenum(ppdeObj)*npt,nt);
for ii=1:size(depvar,2)
    actdepvar=depvar(:,ii);
    tmp=actdepvar(coeffidx).';
    costateval(:,ii)=tmp(:);
end

