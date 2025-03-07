function ocObj=changeparametervalue(ocObj,varargin)
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
% changing one parameter would be ocObj=changeparametervalue(ocObj,1,0.05),
% for changing more parameters 
% ocObj=changeparametervalue(ocObj,[1 2 3],[0.5 1.5 4.02]). 
%
% CHANGEPARAMETERVALUE(OCOBJ,NAME,VALUE) assigns a certain value VALUE to
% the parameter with the name NAME and returns the oc model with this
% changed parameter. One or more parameters can be changed, names of
% different parameters have to be seperated by a comma e.g.
% ocObj=changeparametervalue(ocObj,'alpha,beta',[1.5 2.1]).
% or
% ocObj=changeparametervalue(ocObj,{'alpha','beta'},[1.5 2.1]).


parname=[];
parindex=[];
newparvalue=[];

if nargin>=4 || nargin==1
    ocmaterror('Wrong number of arguments.')  
end
ocObj=cleanmodel(ocObj);
if isnumeric(varargin{1}) && nargin==2      % Second argument is a vector with the values of the new parameters
    ocObj=setparametervalue(ocObj,varargin{1});
    return
end
[oldparvalue totalparname]=parametervalue(ocObj);
if nargin==2 && isocmodel(varargin{1})
    [parvalue parname]=parametervalue(varargin{1});
    newparvalue=oldparvalue;
    for ii=1:length(totalparname)
        cmp=strcmp(totalparname{ii},parname);
        if any(cmp)
            newparvalue(ii)=parvalue(cmp==1);
        end
    end
    ocObj=setparametervalue(ocObj,newparvalue);
    return
end
if nargin>=3
    newparvalue=varargin{2};
end
if ~isnumeric(newparvalue)
    ocmaterror('Third argument must be numeric.')
end
numparvalue=numel(newparvalue);

% Evaluate the input arguments
if nargin>=2
    if ischar(varargin{1})                                      % Second argument is name of a parameter
        varargin{1}=regexprep(varargin{1},'\[||\]','');
        parname=strtrim(regexp(varargin{1},',','split'));
    elseif iscellstr(varargin{1})
        parname=varargin{1};
    elseif nargin==3&&isnumeric(varargin{1})               % Second argument is the index number of a parameter
        parindex=varargin{1};
    else
        ocmaterror('Wrong format of the second argument.')
    end
end

if ~isempty(parname)
    if numel(parname)~=numparvalue
        ocmaterror('Number of parameter arguments and values are different.')
    end
    parindex=zeros(1,numparvalue);
    for ii=1:numparvalue
        idx=find(strcmp(totalparname,parname{ii}));
        if ~isempty(idx)
            parindex(ii)=idx;
        else
            ocmatmsg('Parameter name ''%s'' is unknown.\n',parname{ii})
        end
    end
end
removeidx=find(~parindex);
parindex(removeidx)=[];
newparvalue(removeidx)=[];
if length(parindex)~=length(newparvalue)
    error('Incorrect number of parameter values.')
end
if ~isempty(parindex) && (max(parindex)>numel(oldparvalue) || min(parindex)<1)
    error('Incorrect parameter index.')
end

% change old parameter values to new values
oldparvalue(parindex)=newparvalue;
ocObj=changeparametervalue(ocObj,oldparvalue);