function b=explicitconnectiontime(mmObj)
%
% EXPLICITCONNECTIONTIME determines if the model explicitly depends on the
% connection time for multistage models.
b=[];
if isempty(mmObj)
    return
end
b=explicitconnectiontime(mmObj.stdocmodel.Model);
