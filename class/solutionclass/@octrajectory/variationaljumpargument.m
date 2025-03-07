function out=variationaljumpargument(ocTrj)

out=[];
sinfo=solverinfo(ocTrj);

if isfield(sinfo,'stateconstraint')
    if sinfo.stateconstraint
        if  isfield(sinfo,'vjumpargument')
            out=sinfo.vjumpargument;
        elseif isfield(sinfo,'ventrytimecoordinate')
            out=sinfo.parameters(sinfo.ventrytimecoordinate);
        end
    end
end

if ~isempty(out)
    out=out(:).';
end