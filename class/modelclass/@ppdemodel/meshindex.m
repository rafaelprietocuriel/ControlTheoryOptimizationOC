function coeffidx=meshindex(ppdeObj,solObj,indextype)
%

if nargin==2
    indextype='depvar';
end
coeffidx=[];
if isempty(ppdeObj) || isempty(solObj)
    return
end

if isppdeprimitive(solObj) ||  isppdetrajectory(solObj) || isppdeasymptotic(solObj)
    npt=numpoints(solObj);
end

switch indextype
    case 'depvar'
        n=statenum(ppdeObj);
        idx=1:2*n*npt;
        coeffidx=reshape(idx,npt,[]).';
    case 'state'
        n=statenum(ppdeObj);
        idx=1:n*npt;
        coeffidx=reshape(idx,npt,[]).';
    case 'costate'
        n=statenum(ppdeObj);
        idx=n*npt+1:2*n*npt;
        coeffidx=reshape(idx,npt,[]).';
    case 'control'
        n=controlnum(ppdeObj);
        idx=1:n*npt;
        coeffidx=reshape(idx,npt,[]).';
end

