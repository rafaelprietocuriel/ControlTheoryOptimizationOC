function addmodelpath(mmObj)
%
% ADDMODELPATH add model and data path to MATLAB search path.
%
% ADDMODELPATH(OCOBJ) add the model and data folder for model OCOBJ to the
% MATLAB path.  
%
% ADDMODELPATH(OCOBJ,FORCE) If FORCE is 1 the model and data folder are
% added to the MATLAB path independently if the folder is already on the
% MATLAB path. If FORCE is 0 the user is asked, if the folder is not on the
% MATLAB path

addmodelpath(modelname(mmObj),1);
