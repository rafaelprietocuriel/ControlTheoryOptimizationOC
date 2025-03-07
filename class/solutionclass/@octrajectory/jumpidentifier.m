function out=jumpidentifier(ocTrj)

out=[];
sinfo=solverinfo(ocTrj);

if isfield(sinfo,'stateconstraint') &&  isfield(sinfo,'entryindex') &&  isfield(sinfo,'jumpid')
    out=[sinfo.jumpid(sinfo.entryindex>0); ...
        find(sinfo.entryindex>0)];
end