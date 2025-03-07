function limSet=limitset(mapPrim)

limSet.y=dependentvar(mapPrim);
limSet.arcarg=arcargument(mapPrim);
limSet.arcposition=arcposition(mapPrim);
limSet.modelname=modelname(mapPrim);
limSet.modelparameter=modelparameter(mapPrim);
limSet=occurve(limSet);
