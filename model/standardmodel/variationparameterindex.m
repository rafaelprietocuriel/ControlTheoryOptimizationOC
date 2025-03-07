function out=variationparameterindex(ocStruct)
%
%
out=[];
if isempty(ocStruct)
    return
end
out=retrievemodelinformation(ocStruct,'variationparameterindex');
out=out.value;