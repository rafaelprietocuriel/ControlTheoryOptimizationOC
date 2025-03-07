function out=objectivevalue(statEx)
out=[];
if isempty(statEx)
    return
end

func=str2func([modelname(statEx) 'ObjectiveFunction']);

out=func(variable(statEx),modelparameter(statEx));