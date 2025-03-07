function n=exogenousfunctionnumber(ocObj,varargin)
%
% EXOGENOUSFUNCTIONNUMBER 

n=[];
if isempty(ocObj)
    return
end
info=retrievemultistagemodelinformation(ocObj.stdocmodel.Model,'exogenousfunctionnum');
n=info.value;
