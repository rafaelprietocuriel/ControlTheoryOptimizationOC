function rmmodelpath(ocObj,force)
%
% RMMODELPATH(OCOBJ) removes the model and data path for model OCOBJ from
% the MATLAB path. 
%
% RMMODELPATH(OCOBJ,FORCE) If FORCE is 1 the model and data folder are
% remove from the MATLAB path, otherwise (0) the user is asked if s/he
% wants to continue.

if nargin==1
    force=0;
end

rmmodelpath(modelname(ocObj),force);