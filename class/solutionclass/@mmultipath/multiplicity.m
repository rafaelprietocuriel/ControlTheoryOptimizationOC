function out=multiplicity(ocMultiPath)

out=[];
if isempty(ocMultiPath.userinfo) || ~isfield(ocMultiPath.userinfo,'multiplicity')
    return
end

out=ocMultiPath.userinfo.multiplicity;