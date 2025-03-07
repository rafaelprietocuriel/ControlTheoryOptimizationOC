function b=isautonomous(ocStruct)
%
%
b=[];
if isempty(ocStruct)
    return
end
switch ocStruct.modeltype
    case 'standardmodel'
        out=retrievemodelinformation(ocStruct,'autonomous');
    case 'multistagemodel'
        out=retrievemultistagemodelinformation(ocStruct,'autonomous');
    case 'ppdemodel'
        out=retrieveppdemodelinformation(ocStruct,'autonomous');
    case 'differentialgame'
        out=retrievedifferentialgameinformation(ocStruct,'autonomous');
end
b=out.value==1;