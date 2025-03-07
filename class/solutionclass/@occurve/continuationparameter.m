function val=continuationparameter(ocCuv)

try
    val=ocCuv.userinfo.continuationparameter;
catch
    val=[];
end