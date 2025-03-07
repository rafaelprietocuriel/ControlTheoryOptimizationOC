function J=jacobian(gdynPrim,varargin)
%
% returns the state-costate jacobian of an equilibrium.

J=[];

if isequilibrium(gdynPrim)
    J=linearization(gdynPrim);
    ocObj=loadmodel(gdynPrim);
    arcid=arcargument(gdynPrim);
    coord=implicitcontrolcoordinate(ocObj,arcid);
    if ~isempty(coord)
        par=modelparameter(gdynPrim);
        if isempty(par)
            return
        end
        ocObj=changeparametervalue(ocObj,par);
        statecostatenum=2*statenum(ocObj);
        dUdX=derivativegeneralizedcontrol(gdynPrim);
        J{1}=J{1}(1:statecostatenum,1:statecostatenum)+J{1}(1:statecostatenum,statecostatenum+(1:length(coord)))*dUdX{1};
    end
else
    ocmatmsg('Not yet implemented for limit-cycles.\n')
end
