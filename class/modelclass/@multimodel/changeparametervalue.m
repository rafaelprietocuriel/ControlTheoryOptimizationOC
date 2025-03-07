function mmObj=changeparametervalue(mmObj,varargin)
%
% CHANGEPARAMETERVALUE returns and stdocmodel with new parameters
%
% CHANGEPARAMETERVALUE(OCOBJ,PAR) returns an optimal control model with
% changed parameter values PAR. The number of parameter values PAR have to
% equal the number of parameter values of 'stdocmodel' OCOBJ. If the
% parameter values are changed the 'Result' field of the new model is
% empty. 
% 
% CHANGEPARAMETERVALUE(OCOBJ,IDX,VALUE) assigns the value VALUE to the
% parameter with  index IDX and returns the 'stdocmodel' OCOBJ with this
% changed parameter. One or more parameters can be changed, an example for
% changing one parameter would be mmObj=changeparametervalue(mmObj,1,0.05),
% for changing more parameters 
% mmObj=changeparametervalue(mmObj,[1 2 3],[0.5 1.5 4.02]). 
%
% CHANGEPARAMETERVALUE(OCOBJ,NAME,VALUE) assigns a certain value VALUE to
% the parameter with the name NAME and returns the oc model with this
% changed parameter. One or more parameters can be changed, names of
% different parameters have to be seperated by a comma e.g.
% mmObj=changeparametervalue(mmObj,'alpha,beta',[1.5 2.1]).
% or
% mmObj=changeparametervalue(mmObj,{'alpha','beta'},[1.5 2.1]).


nummod=numberofmodels(mmObj);
mmObj=cleanmodel(mmObj);
for ii=2:nargin
    if ischar(varargin{ii-1}) ||  isnumeric(varargin{ii-1})
        varargin{ii-1}=repmat({varargin{ii-1}},1,nummod);
    end
end
for ii=1:nummod
    switch nargin
        case 2
            mmObj.Model{ii}=changeparametervalue(mmObj.Model{ii},varargin{1}{ii});
        case 3
            mmObj.Model{ii}=changeparametervalue(mmObj.Model{ii},varargin{1}{ii},varargin{2}{ii});
    end
end
