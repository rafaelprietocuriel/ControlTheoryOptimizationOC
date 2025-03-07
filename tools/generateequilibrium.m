function ocEP=generateequilibrium(X,arcid,par,modelname,etype,simpleflag,modeltype)

if nargin<=6
    modeltype='stdocmodel';
end
if nargin<=5
    simpleflag=[];
end
if nargin==4
    simpleflag=[];
end
if isempty(etype)
    etype='dynprimitive';
end
if isempty(simpleflag)
    simpleflag=0;
end
modeltype=str2func(modeltype);
ocObj=modeltype(modelname,[],[],0);
ocObj=changeparametervalue(ocObj,par);

switch etype
    case 'dynprimitive'
        hatx.x=0;
        hatx.y=X;
        hatx.arcarg=arcid;
        hatx.arcinterval=[0 1];
        hatx.modelname=modelname;
        hatx.modelparameter=par;
        if ~simpleflag
            hatx.linearization=jacobian(ocObj,hatx);
        end
        ocEP=dynprimitive(hatx);
    case 'gdynprimitive'
        hatx.odenum=length(X);
        hatx.y{1}=X;
        hatx.x=0;
        hatx.arcarg=arcid;
        hatx.arcinterval=[0 1];
        hatx.modelname=modelname;
        hatx.modelparameter=par;
        ocEP=gdynprimitive(hatx);
        if ~simpleflag
            hatx.linearization=linearization(ocEP);
            ocEP=gdynprimitive(hatx);
        end
    otherwise
        ocEP=[];
        return
end