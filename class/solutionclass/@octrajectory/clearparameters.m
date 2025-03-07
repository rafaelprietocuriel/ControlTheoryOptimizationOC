function ocTrj=clearparameters(ocTrj,paridx)

if nargin==1
    paridx=[];
end
if isempty(paridx) || isempty(parameters(ocTrj))
    ocTrj.solverinfo.parameters=[];
    return
end
if ~all(ismember(paridx,1:length(ocTrj.solverinfo.parameters)))
    ocmaterror('Index exceeds parameter values.')
end
ocTrj.solverinfo.parameters(paridx)=[];