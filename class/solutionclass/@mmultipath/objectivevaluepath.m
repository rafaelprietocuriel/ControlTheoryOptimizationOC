function out=objectivevaluepath(ocMultiPath)

out=[];
if isempty(ocMultiPath.solverinfo)
    return
end

out=cell(1,ocMultiPath.parts);
for ii=1:ocMultiPath.parts
    if isfield(ocMultiPath.solverinfo{ii},'objectivevaluecoord')
        out{ii}=ocMultiPath.solutionclass{ii}.y(ocMultiPath.solverinfo{ii}.objectivevaluecoord,:);
    elseif isfield(ocMultiPath.solverinfo{ii},'objectivevaluecoordinate')
        out{ii}=ocMultiPath.solutionclass{ii}.y(ocMultiPath.solverinfo{ii}.objectivevaluecoordinate,:);
    end
end

if isempty(out{1})
    out=[];
end
