function [val idx]=bifurcationparameter(ocCuv)

try
    val=ocCuv.userinfo.varyparametervalue;
    idx=ocCuv.userinfo.varyparameterindex;
catch
    val=[];
    idx=[];
end