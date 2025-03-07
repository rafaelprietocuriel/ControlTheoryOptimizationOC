function out=lagrangevalue(statEx,varargin)
out=[];
if isempty(statEx)
    return
end

func=str2func([modelname(statEx) 'LagrangeFunction']);

out=func(variable(statEx),lagrangemultiplier(statEx),modelparameter(statEx));