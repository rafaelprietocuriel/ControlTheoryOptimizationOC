function [x,fval,exitflag,info,grad]=ocmatoptimize(ocObj,ctrl,optimStruct,varargin)
%
%
x=[];
fval=[];
exitflag=[];
info=[];
par=parametervalue(ocObj);
gradmethod=getocoptions(optimStruct.option,'GENERAL','GradMethod');
mn=modelname(ocObj);
optimizationfile=str2func([mn optimStruct.optimizationfile]);
gradientfile=str2func([mn optimStruct.gradientfile]);
switch gradmethod
    case 'fmincon'
        [x,fval,exitflag,info,lambda,grad]=fmincon(@minfunction,ctrl,optimStruct.mixedconstraintA,optimStruct.mixedconstraintb,[],[], ...
            optimStruct.lowerbound,optimStruct.upperbound,[],optimStruct.option.STATICOPTIM);
end
    function varargout=minfunction(u)
        X=[optimStruct.statecostate;u];
        varargout{1}=-optimizationfile(optimStruct.time,X,par,optimStruct.arcarg);
        if nargout>1
            varargout{2}=-gradientfile(optimStruct.time,X,par,optimStruct.arcarg).';
        end
    end
end
