function out=constraint(statEx,varargin)
out=[];
if isempty(statEx)
    return
end

func=str2func([modelname(statEx) 'Constraint']);

out=func(variable(statEx),modelparameter(statEx));