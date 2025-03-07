function idx=spatialcontrol(ocStruct)
%
%
idx=[];
if isempty(ocStruct)
    return
end
switch ocStruct.modeltype
    case 'ppdemodel'
        out=retrieveppdemodelinformation(ocStruct,'spatialcontroldependenceindex');
end
idx=out.value;