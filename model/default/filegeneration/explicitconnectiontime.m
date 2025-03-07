function b=explicitconnectiontime(ocStruct)
%
%
b=[];
if isempty(ocStruct)
    return
end
switch ocStruct.modeltype
    case 'standardmodel'
        out=0;
    case 'multistagemodel'
        out=retrievemultistagemodelinformation(ocStruct,'explicitconnectiontime');
end
b=out.value==1;