function out=entryindex(ocTrj)

out=[];
sinfo=solverinfo(ocTrj);

if isfield(sinfo,'stateconstraint') &&  isfield(sinfo,'entryindex')
    out=[find(sinfo.entryindex>0); ...
        sinfo.entryindex(find(sinfo.entryindex>0))];
end