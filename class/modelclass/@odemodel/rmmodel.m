function varargout=rmmodel(odeObj,force)
%
% RMMODEL removes a model from OCMat
%
% RMMODEL(OCOBJ) deletes the model files of model OCOBJ and removes
% the model folders. Additionally the folders are removed from the MATLAB
% path
%
% RMMODEL(OCOBJ,FORCE) default value FORCE=0
% FORCE=1: forces to remove the folders
% without asking.
% FORCE=0: the user is asked if s/he wants to continue.
%
% [STAT,MESS,ID]=RMMODEL(...) Here, STAT is 1 for success and is 0 for
% error, and message MESS, messageid ID.

if nargin==1
    force=0;
end
[varargout{1:nargout}]=rmmodel(modelname(odeObj),force);