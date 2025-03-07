function out=partstructure(ocMultiPath)

out=[];
if isempty(ocMultiPath.userinfo) || ~isfield(ocMultiPath.userinfo,'partstructure')
    return
end

out=ocMultiPath.userinfo.partstructure;