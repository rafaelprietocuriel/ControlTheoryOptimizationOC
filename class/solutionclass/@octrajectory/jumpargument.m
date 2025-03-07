function out=jumpargument(ocTrj)

out=[];
sinfo=solverinfo(ocTrj);

if isfield(sinfo,'stateconstraint')
    if sinfo.stateconstraint
        if isfield(sinfo,'entrytimecoordinate')
            out=sinfo.parameters(sinfo.entrytimecoordinate);
        elseif isfield(sinfo,'jumpargument')
            out=sinfo.jumpargument;
        end
    end
end

if ~isempty(out)
    out=out(:).';
end