function out=variationalfreetime(ocTrj)

out=[];
sinfo=solverinfo(ocTrj);

if isfield(sinfo,'vfreetimecoord')
    out=sinfo.parameters(sinfo.vfreetimecoord);
elseif isfield(sinfo,'vfreetime')
    out=sinfo.vfreetime;
end

if ~isempty(out)
    out=out(:).';
end