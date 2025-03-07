function b=variationparameter(ocStruct)
%
%
b=[];
if isempty(ocStruct)
    return
end
out=retrievemodelinformation(ocStruct,'variationparameternum');
b=out.value~=0;