function ocStruct=odemodelsupplement(ocStruct,supplementtype)
%
%
switch supplementtype
    case 'maininitfile'
        % if no arc is defined set unconstrained case
        if ~isfield(ocStruct,'arc') || ~isfield(ocStruct.arc,'identifier') ...
                || isempty(ocStruct.arc.identifier)
            ocStruct.arc.identifier{1}='0';
            ocStruct.arc.argument=0;
            ocStruct.arc.num=1;
        end

        % set default variable name for independent variable
        if ~isfield(ocStruct.variable,'independent') || ~isfield(ocStruct.variable.independent,'name') ...
                || isempty(ocStruct.variable.independent.name)
            ocStruct.variable.independent.name=getbasicname('independent');
            ocStruct.variable.independent.num=1;
        end

        % determine if the problem is autonomous, neither the objective function
        % nor the state dynamics depend explicitly on the independent variable
        ocStruct.variable.independent.property.autonomous=testautonomous(ocStruct);
end

function autflag=testautonomous(ocStruct)
%
% the OC problem is called autonomous if the integrand of the objective
% function, state dynamics and constraint functions do not explicitly
% depend on time. Furthermore it is assumed that the discount factor is
% given by the exponential function. Under these assumptions the canonical
% system can be formulated as an autonomous system of ODEs justifying the
% terminology
autflag=1;
if isempty(ocStruct)
    return
end
independentvar=ocStruct.variable.independent.name;
for ii=1:numel(ocStruct.constraint.derivative.state.term)
    if ~isempty(ocStruct.constraint.derivative.state.term{ii})
        autflag=autflag && isempty(regexp(ocStruct.constraint.derivative.state.term{ii},['\<' independentvar '\>']));
    end
end
if isfield(ocStruct,'exogenousfunction')
    fn=fieldnames(ocStruct.exogenousfunction);
    for ii=1:length(fn)
            autflag=autflag && isempty(regexp(ocStruct.exogenousfunction.(fn{ii}).term,['\<' independentvar '\>']));
    end
end