function idx=nonspatialcontrol(ocStruct)
%
%
idx=[];
if isempty(ocStruct)
    return
end
switch ocStruct.modeltype
    case 'ppdemodel'
        out=retrieveppdemodelinformation(ocStruct,'nonspatialcontroldependenceindex');
end
idx=out.value;