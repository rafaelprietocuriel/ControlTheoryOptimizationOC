function limSet=limitset(dynPrim)

limSet.y=dependentvar(dynPrim);
limSet.arcarg=arcargument(dynPrim);
limSet.arcposition=arcposition(dynPrim);
limSet.modelname=modelname(dynPrim);
limSet.modelparameter=modelparameter(dynPrim);
limSet=occurve(limSet);
