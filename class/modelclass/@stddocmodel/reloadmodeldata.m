function varargout=reloadmodeldata(ocObj)
%
% LOAD loads the data for oc model into the workspace
%
% LOAD(OCOBJ) loads the data for oc model into the workspace
%
% LOAD(OCOBJ,'WORKSPACE') loads data previously stored by the
% SAVE(OCOBJ,'WORKSPACE') command.
%
% LOAD(OCOBJ,'MODELDATA') loads model structure of the model into the
% workspace and variable name 'ocStruct'. 
%
% LOAD(OCOBJ,FORMAT) uses the format string FORMAT (see SPRINTF for
% details) for the parameter values used for the file name generation.
%
% LOAD(OCOBJ,FORMAT,FORCE) FORCE = 1 forces to overwrite the results of
% OCOBJ by the data stored from a previous session. FORCE = 0 if a data
% file exists the user is asked if s/he wants to proceed and overwrite the
% model data.
%
% LOAD(OCOBJ,FORMAT,FORCE,FN) FN provides an alternative filename.
%
% LOAD(OCOBJ,FORMAT,FORCE,IDX) filename is generated of the model
% name and the parameter values for the IDX'th parameter values. 
%
% LOAD(OCOBJ,FORMAT,FORCE,FN,SAVEDIR) SAVEDIR provides an alternative
% folder. 

ocObj.Model=loadmodeldata(modelname(ocObj));
assignin('caller',inputname(1),ocObj);
if nargout==1
    varargout{1}=true;
end
