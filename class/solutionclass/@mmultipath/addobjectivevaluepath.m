function ocMultiPath=addobjectivevaluepath(ocMultiPath,objpath,objcoord)

if isempty(ocMultiPath.userinfo) || ~isfield(ocMultiPath.userinfo,'partstructure')
    return
end

for ii=1:ocMultiPath.parts
    ocMultiPath.solutionclass{ii}=addobjectivevaluepath(ocMultiPath.solutionclass{ii},objpath{ii},objcoord);
    ocMultiPath.solverinfo{ii}=solverinfo(ocMultiPath.solutionclass{ii});
end