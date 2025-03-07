function varargout=reloadmodeldata(ocObj)
%
% RELOADMODELDATA reloads the model data structure after a possible change
% of the model formulation. The original parameter values are kept.
%
% RELOADMODELDATA(OCOBJ) reloads the model data structure into the model of
% the calling workspace.

if ~isempty(ocObj)
    ocObjTmp=stdocmodel(modelname(ocObj));
    ocObjTmp=changeparametervalue(ocObjTmp,ocObj);
    oldparvalue=parametervalue(ocObjTmp);
    ocObj.Model=loadmodeldata(modelname(ocObj));
end
ocObj=setparametervalue(ocObj,oldparvalue);
assignin('caller',inputname(1),ocObj);
if nargout==1
    varargout{1}=true;
end
