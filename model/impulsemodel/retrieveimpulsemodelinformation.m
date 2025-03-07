function data=retrieveimpulsemodelinformation(ocStruct,propertyclass,arcidentifier,symkernel)
%
% RETRIEVEDIFFMODELINFORMATION returns a structure with information about the
% ocmat model
% this function is the interface to the structure OCSTRUCT as derived
% from the models initialization file. Its purpose is to allow the same
% commands even if the structure OCSTRUCT is changed. In that case only
% 'retrieveodemodelinformation' has to be adapted.

data.type='';
data.value='';
data.description='';

if isempty(ocStruct)
    return
end
if nargin==2
    arcidentifier='';
end

switch propertyclass
    case 'arcidentifier'
        % arc identifiers specify different combinations of (in)active
        % constraints and/or multiple solutions of the  Hamiltonian
        % maximizing condition
        data.type='char';
        data.value=ocStruct.arc.identifier;
        data.description='identifier to differentiate between specifications of the canonical system';

    case 'argument'
        % arcidentifiers (char) are transformed the corresponding arguments
        % (numeric)
        data.type='integervector';
        data.value=ocStruct.arc.argument;
        data.description='arcarguments is a vector containing the numeric values of the arcidentifiers';

    case 'constraintcombination'
        % allowed constraint combinations specified by the user in the
        % initialization file
        arcarg=arcidentifier2arcindex(arcidentifier);
        data.type='cellchar';
        data.value=regexp(strtrim(ocStruct.arc.constraintcombination{arcarg}),'\s*','split');
        data.description='allowed constraint combinations';

    case 'constraintcombinationindex'
        % zero one vector representing constraint combination
        data.type='integervector';
        constraintcombination=retrieveimpulsemodelinformation(ocStruct,'constraintcombination',arcidentifier);
        inequalitycontrolconstraintidentifier=retrieveimpulsemodelinformation(ocStruct,'inequalitycontrolconstraintidentifier');
        inequalitycontrolconstraintnum=retrieveimpulsemodelinformation(ocStruct,'inequalitycontrolconstraintnum');
        idx=zeros(1,inequalitycontrolconstraintnum.value);
        for ii=1:numel(constraintcombination.value)
            idx(strcmp(inequalitycontrolconstraintidentifier.value,constraintcombination.value{ii}))=1;
        end
        data.value=idx;
        data.description='zero one vector representing constraint combination';

    case 'independent'
        % variable name of the indepent variable (usually time)
        data.type='char';
        data.value=ocStruct.variable.independent.name;
        data.description='independent variable';

%     case 'impulsetime'
%         % variable name of the indepent variable (usually time)
%         data.type='char';
%         data.value=ocStruct.variable.impulsetime.name;
%         data.description='impulsetime variable';

    case 'impulsetime'
        % variable name of the indepent variable (usually time)
        data.type='char';
        data.value=ocStruct.variable.impulsetime.name;
        data.description='independent variable';

    case 'endtime'
        % variable name of the indepent variable (usually time)
        data.type='char';
        data.value=ocStruct.variable.endtime.name;
        data.description='endtime variable';

    case 'autonomous'
        % variable name of the indepent variable (usually time)
        data.type='boolean';
        data.value=ocStruct.variable.independent.property.autonomous;
        data.description='model not explicitly depending on time';

    case 'statename'
        % variable name of the state(s) as defined by the user in the
        % initialization file
        data.type='cellchar';
        data.value=ocStruct.variable.state.name;
        data.description='name of state';

    case 'statenamel'
        % variable name of the state(s) as defined by the user in the
        % initialization file
        data.type='cellchar';
        for ii=1:length(ocStruct.variable.state.name)
            data.value{ii}=[ocStruct.variable.state.name{ii} 'L'];
        end
        data.description='name of state';

    case 'statenamer'
        % variable name of the state(s) as defined by the user in the
        % initialization file
        data.type='cellchar';
        for ii=1:length(ocStruct.variable.state.name)
            data.value{ii}=[ocStruct.variable.state.name{ii} 'R'];
        end
        data.description='name of state';

    case 'costatename'
        % variable name of the costate(s) as defined by the user in the
        % initialization file
        data.type='cellchar';
        data.value=ocStruct.variable.costate.name;
        data.description='name of costate';

    case 'costatenamel'
        % variable name of the costate(s) as defined by the user in the
        % initialization file
        data.type='cellchar';
        data.value=cell(1:length(ocStruct.variable.costate.name));
        for ii=1:length(ocStruct.variable.costate.name)
            data.value{ii}=[ocStruct.variable.costate.name{ii} 'L'];
        end
        data.description='name of costate';

    case 'costatenamer'
        % variable name of the costate(s) as defined by the user in the
        % initialization file
        data.type='cellchar';
        data.value=cell(1:length(ocStruct.variable.costate.name));
        for ii=1:length(ocStruct.variable.costate.name)
            data.value{ii}=[ocStruct.variable.costate.name{ii} 'R'];
        end
        data.description='name of costate';

    case 'icontrolname'
        % variable name of the control(s) as defined by the user in the
        % initialization file
        data.type='cellchar';
        data.value=ocStruct.variable.icontrol.name;
        data.description='name of impulse control';

    case 'controlname'
        % variable name of the control(s) as defined by the user in the
        % initialization file
        data.type='cellchar';
        try
            data.value=ocStruct.variable.control.name;
        catch
            data.value='';
        end
        data.description='name of control';


    case 'exogenousfunctionnum'
        % variable name of the control(s) as defined by the user in the
        % initialization file
        data.type='cellchar';
        if isfield(ocStruct,'exogenousfunction')
            data.value=length(fieldnames(ocStruct.exogenousfunction));
        else
            data.value=0;
        end
        data.description='names of exogenous functions';

    case 'exogenousfunctionname'
        % variable name of the control(s) as defined by the user in the
        % initialization file
        data.type='cellchar';
        if isfield(ocStruct,'exogenousfunction')
            data.value=fieldnames(ocStruct.exogenousfunction).';
        else
            data.value=[];
        end
        data.description='names of exogenous functions';

        
    case 'exogenousfunctionnamewithargument'
        % variable name of the control(s) as defined by the user in the
        % initialization file
        data.type='cellchar';
        if isfield(ocStruct,'exogenousfunction')
            exogenousfunctionname=retrieveimpulsemodelinformation(ocStruct,'exogenousfunctionname');
            exogenousfunctionargument=retrieveimpulsemodelinformation(ocStruct,'exogenousfunctionargument');
             for ii=1:length(exogenousfunctionname.value)
                data.value{ii}=[exogenousfunctionname.value{ii} '(' exogenousfunctionargument.value{ii} ')'];
            end
       else
            data.value=[];
        end
        data.description='names of exogenous functions with arguments';

    case 'exogenousfunctionargument'
        % variable name of the control(s) as defined by the user in the
        % initialization file
        data.type='cellchar';
        exogenousfunctionname=retrieveimpulsemodelinformation(ocStruct,'exogenousfunctionname');
        if ~isempty(exogenousfunctionname.value)
            for ii=1:length(exogenousfunctionname.value)
                data.value{ii}=ocStruct.exogenousfunction.(exogenousfunctionname.value{ii}).argument;
            end
        else
            data.value=[];
        end
        data.description='arguments of exogenous functions';

    case 'exogenousfunctionnamedx'
        % variable name of the control(s) as defined by the user in the
        % initialization file
        data.type='cellmatrixchar';
        exogenousfunctionname=retrieveimpulsemodelinformation(ocStruct,'exogenousfunctionname');
        exogenousfunctionargument=retrieveimpulsemodelinformation(ocStruct,'exogenousfunctionargument');
        statename=retrieveimpulsemodelinformation(ocStruct,'statename');
        statenum=retrieveimpulsemodelinformation(ocStruct,'statenum');
        for ii=1:length(exogenousfunctionname.value)
            if ~isempty(exogenousfunctionargument.value{ii})
                counter=0;
                while counter<statenum.value
                    counter=counter+1;
                    if ~isempty(regexp(ocStruct.exogenousfunction.(exogenousfunctionname.value{ii}).argument,['\<' statename.value{counter} '\>'],'once'))
                        break
                    end
                end
                if counter<=statenum.value
                    data.value{ii}=['D1' exogenousfunctionname.value{ii} '_D' exogenousfunctionargument.value{ii} '(' exogenousfunctionargument.value{ii} ')'];
                else
                    data.value{ii}='';
                end
                
            else
                data.value='';
            end
        end
        data.description='arguments of exogenous functions';

    case 'exogenousfunctionnamederivatived2x2'
        % variable name of the control(s) as defined by the user in the
        % initialization file
        data.type='cellchar';
        exogenousfunctionname=retrieveimpulsemodelinformation(ocStruct,'exogenousfunctionname');
        exogenousfunctionargument=retrieveimpulsemodelinformation(ocStruct,'exogenousfunctionargument');
        statename=retrieveimpulsemodelinformation(ocStruct,'statename');
        statenum=retrieveimpulsemodelinformation(ocStruct,'statenum');
        for ii=1:length(exogenousfunctionname.value)
            if ~isempty(exogenousfunctionargument.value{ii})
                counter=0;
                while counter<statenum.value
                    counter=counter+1;
                    if ~isempty(regexp(ocStruct.exogenousfunction.(exogenousfunctionname.value{ii}).argument,['\<' statename.value{counter} '\>'],'once'))
                        break
                    end
                end
                if counter<=statenum.value
                    data.value{ii}=['D2' exogenousfunctionname.value{ii} 'Dx2'];
                else
                    data.value{ii}='';
                end
                
            else
                data.value='';
            end
        end
        data.description='arguments of exogenous functions';

    case 'exogenousfunctionterm'
        % variable name of the control(s) as defined by the user in the
        % initialization file
        data.type='cellchar';
        exogenousfunctionname=retrieveimpulsemodelinformation(ocStruct,'exogenousfunctionname');
        if ~isempty(exogenousfunctionname.value)
            for ii=1:length(exogenousfunctionname.value)
                data.value{ii}=ocStruct.exogenousfunction.(exogenousfunctionname.value{ii}).term;
            end
        else
            data.value=[];
        end
        data.description='terms of exogenous functions';

    case 'parametername'
        % variable name(s) of the parameter(s)
        data.type='cellchar';
        data.value=fieldnames(ocStruct.parameter.variable);
        data.description='name of parameter variable';

    case 'parametervalue'
        % vector of user provided parameter values
        data.type='double';
        parametername=retrieveimpulsemodelinformation(ocStruct,'parametername');
        parametervaluetermindex=retrieveimpulsemodelinformation(ocStruct,'parametervaluetermindex');
        parameteralgebraictermindex=retrieveimpulsemodelinformation(ocStruct,'parameteralgebraictermindex');
        data.value=zeros(1,numel(parametername.value));
        for ii=parametervaluetermindex.value
            data.value(ii)=ocStruct.parameter.variable.(parametername.value{ii});
            if ~isempty(parameteralgebraictermindex.value)
                eval([parametername.value{ii} '=ocStruct.parameter.variable.' parametername.value{ii} ';'])
            end
        end
        %             term=ocStruct.parameter.variable.(parametername.value{ii});
        for ii=parameteralgebraictermindex.value
            data.value(ii)=eval(ocStruct.parameter.variable.(parametername.value{ii}));
            eval([parametername.value{ii} '=data.value(ii);'])
        end
        %             data.value(ii)=double(term);
        data.description='vector of user provided parameter values';

    case 'discountrate'
        data.type='double';
        parametername=retrieveimpulsemodelinformation(ocStruct,'parametername');
        parametervalue=retrieveimpulsemodelinformation(ocStruct,'parametervalue');
        discountratevariable=retrieveimpulsemodelinformation(ocStruct,'discountratevariable');
        data.value=parametervalue.value(strcmp(parametername.value,discountratevariable.value));
        data.description='value of the discount rate';

    case 'variablename'
        % general variable name(s) of the optimal control problem
        data.type='cellchar';
        data.value=fieldnames(ocStruct.variable);
        data.description='general variable name';

    case 'lagrangemultipliercontrolname'
        % variable name of the Lagrange multiplier(s) for inequality
        % (control) constraints as defined by the user in the
        % initialization file
        if isfield(ocStruct.variable,'lagrangemultcc')
            data.type='cellchar';
            data.value=ocStruct.variable.lagrangemultcc.name;
        end
        data.description='name of Lagrange multiplier';

    case 'inequalitycontrolconstraintidentifier'
        if isfield(ocStruct.constraint,'function') && ...
                isfield(ocStruct.constraint.function,'control') && ...
                ocStruct.constraint.function.control.num
            data.value=ocStruct.constraint.function.control.identifier;
            data.type='cellchar';
        end

    case 'statenum'
        % the number of state(s)
        data.type='integer';
        data.value=ocStruct.variable.state.num;
        data.description='number of state';

    case 'costatenum'
        % the number of costate(s)
        data.type='integer';
        data.value=ocStruct.variable.costate.num;
        data.description='number of costate';

    case 'icontrolnum'
        % the number of control(s)
        data.type='integer';
        data.value=ocStruct.variable.icontrol.num;
        data.description='number of impulse control';

    case 'controlnum'
        % the number of control(s)
        data.type='integer';
        if isfield(ocStruct.variable,'control')
            data.value=ocStruct.variable.control.num;
        else
            data.value=00;
        end
        data.description='number of control';

    case 'arcnum'
        % the number of different arcs, corresponding to different
        % functional specifications of the canonical system
        data.type='integer';
        data.value=ocStruct.arc.num;
        data.description='number of different arcs';
    case 'jumpidnum'
        % the number of different arcs, corresponding to different
        % functional specifications of the canonical system
        data.type='integer';
        data.value=0;%ocStruct.arc.num;
        data.description='number of different arcs';

    case 'parameternum'
        % number of parameter(s)
        data.type='integer';
        if isfield(ocStruct,'parameter')
            data.value=ocStruct.parameter.num;
        else
            data.value=0;
        end
        data.description='number of exogenous parameters';

    case 'parameteralgebraictermindex'
        % number of parameter(s)
        data.type='integer vector';
        if isfield(ocStruct,'parameter')
            if isfield(ocStruct.parameter,'algebraictermidx') && ~isempty(ocStruct.parameter.algebraictermidx)
                data.value=ocStruct.parameter.algebraictermidx;
            else
                data.value=[];
            end
        else
            data.value=[];
        end
        data.description='index of parameters given by an algebraic equation';

    case 'parametervaluetermindex'
        % number of parameter(s)
        data.type='integer vector';
        parameternum=retrieveimpulsemodelinformation(ocStruct,'parameternum');
        parameteralgebraictermindex=retrieveimpulsemodelinformation(ocStruct,'parameteralgebraictermindex');
        if isfield(ocStruct,'parameter')
            data.value=setdiff(1:parameternum.value,parameteralgebraictermindex.value);
        else
            data.value=[];
        end
        data.description='index of parameters given by an explicit value';

    case 'parametervalueexpression'
        % vector of user provided parameter values
        data.type='cellchar';
        parametername=retrieveimpulsemodelinformation(ocStruct,'parametername');
        for ii=1:numel(parametername.value)
            if ~isnumeric(ocStruct.parameter.variable.(parametername.value{ii}))
                data.value{ii}=ocStruct.parameter.variable.(parametername.value{ii});
            else
                data.value{ii}=[];
            end
        end
        data.description='vector of user provided parameter value terms';
    case 'odedim'
        % the number of state(s)
        data.type='integer';
        data.value=ocStruct.variable.state.num+ocStruct.variable.costate.num;
        data.description='number of ODEs for the canonical system';

    case 'inequalitycontrolconstraint'
        % inequality constraints which explicitly include control variables
        % (and state variables).
        if isfield(ocStruct.constraint,'function') &&  ...
                isfield(ocStruct.constraint.function,'control') && ...
                ocStruct.constraint.function.control.num

            data.value=ocStruct.constraint.function.control.term;
            data.type='mathchar';
        end
        data.description='mixed inequality constraints';

    case 'inequalitycontrolconstraintnum'
        % number of inequality constraints explicitly including control
        % variables (and state variables).
        data.type='integer';
        if isfield(ocStruct.constraint,'function') && isfield(ocStruct.constraint.function,'control')
            data.value=ocStruct.constraint.function.control.num;
        else
            data.value=0;
        end
        data.description='number of mixed inequality constraints';


    case 'inequalitystateconstraint'
        % inequality constraints which only include state variables
        if isfield(ocStruct.constraint,'function') &&  ...
                isfield(ocStruct.constraint.function,'state') && ...
                ocStruct.constraint.function.state.num

            data.value=ocStruct.constraint.function.state.term;
            data.type='mathchar';
        end
        data.description='pure state inequality constraints';

    case 'inequalitystateconstraintnum'
        % number of inequality constraints explicitly including only state
        % variables.
        data.type='integer';
        if isfield(ocStruct.constraint,'function') && isfield(ocStruct.constraint.function,'state')
            data.value=ocStruct.constraint.function.state.num;
        else
            data.value=0;
        end
        data.description='number of pure state inequality constraints';

    case 'algebraicequationnum'
        maximizingimplicitderivativevariable=retrieveimpulsemodelinformation(ocStruct,'maximizingimplicitderivativevariable',arcidentifier);
        data.type='integer';
        data.value=0;%numel(maximizingimplicitderivativevariable.value);
        data.description='number of algebraic equations';

    case 'totalalgebraicequationnum'
        totalimplicitvariable=retrieveimpulsemodelinformation(ocStruct,'totalimplicitvariable');
        data.type='integer';
        data.value=numel(totalimplicitvariable.value);
        data.description='number of algebraic equations';

    case 'canonicalsystemequationnum'
        statenum=retrieveimpulsemodelinformation(ocStruct,'statenum');
        costatenum=retrieveimpulsemodelinformation(ocStruct,'costatenum');
        data.type='integer';
        data.value=statenum.value+costatenum.value;
        data.description='total number of algebraic differential equations';

    case 'statedynamics'
        % ODEs describing the evolution of the state(s)
        data.type='mathchar';
        data.value=ocStruct.constraint.derivative.state.term;
        data.description='ODEs of the state dynamics';

    case 'stateevent'
        % ODEs describing the evolution of the state(s)
        data.type='mathchar';
        data.value=ocStruct.constraint.event.state.term;
        data.description='Events of the state dynamics';

    case 'objectivetype'
        % type of the terms which build up the objective value function
        % (usually sum and Salvage value)
        data.type='mathchar';
        data.value=fieldnames(ocStruct.objective);
        data.description='type of terms which build up the objective value function';

    case 'objectiveintegrand'
        % summand for objective value function without the possible
        % discounting factor
        data.type='mathchar';
        data.value=ocStruct.objective.integral.function.term;
        data.description='summand for objective value function';

    case 'objectivesummand'
        % summand for objective value function without the possible
        % discounting factor
        data.type='mathchar';
        data.value=ocStruct.objective.sum.function.term;
        data.description='summand for objective value function';

    case 'discountedobjectivesummand'
        % summand for objective value function without the possible
        % discounting factor
        data.type='mathchar';
        data.value=[ocStruct.objective.sum.discountfactor.term '*(' ocStruct.objective.sum.function.term ')'];
        data.description='summand for objective value function';

    case 'salvagevalue'
        % Salvagevalue
        data.type='mathchar';
        if isfield(ocStruct.objective,'sum') && isfield(ocStruct.objective.sum,'endtime')
            data.value=ocStruct.objective.sum.endtime.function.term;
            data.description='salvagevalue';
        end

    case 'discountratevariable'
        % discount variable used in the discount factor of the objective
        % value function
        if isfield(ocStruct.objective,'integral')
            data.type='char';
            data.value=ocStruct.objective.integral.discountrate;
            data.description='discount variable of the objective value function';
        end

    case 'discountfactor'
        % discount variable used in the discount factor of the objective
        % value function
        if isfield(ocStruct.objective,'integral')
            data.type='mathchar';
            data.value=ocStruct.objective.integral.discountfactor.term;
            data.description='discount factor of the objective value function';
        end

    case 'endtimediscountratevariable'
        % discount variable used in the discount factor of the objective
        % value function
        if isfield(ocStruct.objective,'sum') && isfield(ocStruct.objective.sum,'endtime')
            data.type='char';
            data.value=ocStruct.objective.sum.endtime.discountrate;
            data.description='discount variable of salvagevalue';
        end

    case 'endtimediscountfactor'
        % discount variable used in the discount factor of the objective
        % value function
        if isfield(ocStruct.objective,'sum')
            data.type='mathchar';
            data.value=ocStruct.objective.sum.endtime.discountfactor.term;
            data.description='discount factor of the salvage value';
        end

    case 'sumdiscountratevariable'
        % discount variable used in the discount factor of the objective
        % value function
        if isfield(ocStruct.objective,'sum') && isfield(ocStruct.objective.sum,'function')
            data.type='char';
            data.value=ocStruct.objective.sum.discountrate;
            data.description='discount variable of objecitve summand';
        end

    case 'sumdiscountfactor'
        % discount variable used in the discount factor of the objective
        % value function
        if isfield(ocStruct.objective,'sum') && isfield(ocStruct.objective.sum,'function')
            data.type='mathchar';
            data.value=ocStruct.objective.sum.discountfactor.term;
            data.description='discount factor of objecitve summand';
        end

    case 'solutionindex'
        % if more than one solution of Hamiltonian maximizing condition
        % exists this index determines which of the solutions sould be kept
        % by default. During the initialization porcess the user is
        % explicitly asked to confirm this selection
        uniqueconstraintcombination=unique(ocStruct.arc.constraintcombination);
        controlvaluecombination=ocStruct.arc.controlvaluecombination;
        if numel(uniqueconstraintcombination)==numel(ocStruct.arc.constraintcombination)
            data.value=controlvaluecombination;
        else
            for ii=1:numel(uniqueconstraintcombination)
                data.value{ii}=controlvaluecombination(strcmp(ocStruct.arc.constraintcombination,uniqueconstraintcombination{ii}));
            end
        end
        data.description='used index for multiple Hamilton maximizing condition solutions';

    case 'discobjectivesummand'
        objectivefunction=retrieveimpulsemodelinformation(ocStruct,'objectivesummand');
        discountratevariable=retrieveimpulsemodelinformation(ocStruct,'discountratevariable');
        independentvariable=retrieveimpulsemodelinformation(ocStruct,'independent');
        if isfield(ocStruct.objective,'sum')
            data.type='mathchar';
            data.value=['exp(-' discountratevariable.value '*' independentvariable.value ')*(' objectivefunction.value ')'];
            data.description='the discounted objective summand';
        end

    case 'discobjectivefunction'
        objectivefunction=retrieveimpulsemodelinformation(ocStruct,'objectiveintegrand');
        discountratevariable=retrieveimpulsemodelinformation(ocStruct,'discountratevariable');
        independentvariable=retrieveimpulsemodelinformation(ocStruct,'independent');
        if isfield(ocStruct.objective,'integral')
            data.type='mathchar';
            data.value=['exp(-' discountratevariable.value '*' independentvariable.value ')*(' objectivefunction.value ')'];
            data.description='the discounted objective integrand';
        end

    case 'discobjectivefunctionderivativetime'
        specificdiscobjectivefunction=retrieveimpulsemodelinformation(ocStruct,'specificdiscobjectivefunction',arcidentifier,symkernel);
        timevariable=retrieveimpulsemodelinformation(ocStruct,'independent');
        variable=cell2vectorstring(timevariable.value);
        data.value=string2cell(removematrixstring(ocmatdiff(specificdiscobjectivefunction.value,variable,symkernel)));
        data.type='mathchar';
        data.description='the discounted objective integrand';

    case 'specificdiscobjectivefunction'
        discobjectivefunction=retrieveimpulsemodelinformation(ocStruct,'discobjectivefunction');
        optimalvalue=retrieveimpulsemodelinformation(ocStruct,'optimalvalue',arcidentifier);
        maximizingvariable=retrieveimpulsemodelinformation(ocStruct,'maximizingvariable');
        for ii=1:numel(maximizingvariable.value)
            discobjectivefunction.value=ocmatsubs(discobjectivefunction.value,[maximizingvariable.value{ii} '=' optimalvalue.value.(maximizingvariable.value{ii})],symkernel);
        end
        data.type='mathchar';
        data.value=discobjectivefunction.value;
        data.description='the discounted objective integrand';

    case 'discobjectivefunctionDX'
        costatename=retrieveimpulsemodelinformation(ocStruct,'costatename');
        statename=retrieveimpulsemodelinformation(ocStruct,'statename');
        statevector=cell2vectorstring([statename.value costatename.value]);
        specificdiscobjectivefunction=retrieveimpulsemodelinformation(ocStruct,'specificdiscobjectivefunction',arcidentifier,symkernel);
        data.type='mathcellchar';
        data.value=string2cell(ocmatjacobian(['[' specificdiscobjectivefunction.value ']'],statevector,symkernel),'vector');
        data.description='derivative of the discounted objective integrand with respect to the state and costate';

    case 'discobjectivefunctionDP'
        specificdiscobjectivefunction=retrieveimpulsemodelinformation(ocStruct,'specificdiscobjectivefunction',arcidentifier,symkernel);
        parametername=retrieveimpulsemodelinformation(ocStruct,'parametername',arcidentifier);
        variable=cell2vectorstring(parametername.value);
        data.type='mathcellchar';
        data.value=string2cell(ocmatjacobian(['[' specificdiscobjectivefunction.value ']'],variable,symkernel),'vector');
        data.description='derivative of the discounted objective integrand with respect to the state and costate';

    case 'discimpulseobjectivefunction'
        objectivefunction=retrieveimpulsemodelinformation(ocStruct,'objectivesummand');
        discountratevariable=retrieveimpulsemodelinformation(ocStruct,'discountratevariable');
        independentvariable=retrieveimpulsemodelinformation(ocStruct,'impulsetime');
        if isfield(ocStruct.objective,'integral')
            data.type='mathchar';
            data.value=['exp(-' discountratevariable.value '*' independentvariable.value ')*(' objectivefunction.value ')'];
            data.description='the discounted objective integrand';
        end

    case 'arcwitheequalpontryaginfunction'
        % in cases where the Hamiltonian maximizing condition yields
        % multiple result different solutions are identified by different
        % arc identifiers. 'arcwitheequalpontryaginfunction' returns the
        % arc idnetifiers which correspond to the same specification of
        % the Pontryagin function (generalized Hamilton function including
        % terms for (inequality) constraints)
        uniqueconstraintcombination=unique(ocStruct.arc.constraintcombination);
        arcidentifier=ocStruct.arc.identifier;
        if numel(uniqueconstraintcombination)==numel(ocStruct.arc.constraintcombination)
            for ii=1:numel(uniqueconstraintcombination)
                data.value{ii}=arcidentifier(ii);
            end
        else
            for ii=1:numel(uniqueconstraintcombination)
                data.value{ii}=arcidentifier(strcmp(ocStruct.arc.constraintcombination,uniqueconstraintcombination{ii}));
            end
        end
        data.type='cellchar';
        data.description='arcidentifer corresponding to same Pontryaginfunction';

    case 'arcwitheequalimpulsepontryaginfunction'
        % in cases where the Hamiltonian maximizing condition yields
        % multiple result different solutions are identified by different
        % arc identifiers. 'arcwitheequalpontryaginfunction' returns the
        % arc idnetifiers which correspond to the same specification of
        % the Pontryagin function (generalized Hamilton function including
        % terms for (inequality) constraints)
        uniqueconstraintcombination=unique(ocStruct.arc.constraintcombination);
        arcidentifier=ocStruct.arc.identifier;
        if numel(uniqueconstraintcombination)==numel(ocStruct.arc.constraintcombination)
            for ii=1:numel(uniqueconstraintcombination)
                data.value{ii}=arcidentifier(ii);
            end
        else
            for ii=1:numel(uniqueconstraintcombination)
                data.value{ii}=arcidentifier(strcmp(ocStruct.arc.constraintcombination,uniqueconstraintcombination{ii}));
            end
        end
        data.type='cellchar';
        data.description='arcidentifer corresponding to same Pontryaginfunction';

    case 'equationvariablename'
        % dependent variables of the canonical system for a specific arc
        data.type='cellchar';
        statename=retrieveimpulsemodelinformation(ocStruct,'statename');
        costatename=retrieveimpulsemodelinformation(ocStruct,'costatename');
        implicitcontrolname=retrieveimpulsemodelinformation(ocStruct,'implicitnonlinearcontrol',arcidentifier);
        data.value=[statename.value costatename.value implicitcontrolname.value];
        data.description='dependent variables of the canonical system';

    case 'equationvariablenamet'
        % dependent variables of the canonical system for a specific arc
        data.type='cellchar';
        statenamet=retrieveimpulsemodelinformation(ocStruct,'statenamet');
        costatenamet=retrieveimpulsemodelinformation(ocStruct,'costatenamet');
        implicitnonlinearcontroltp1=retrieveimpulsemodelinformation(ocStruct,'implicitnonlinearcontroltp1',arcidentifier);
        data.value=[statenamet.value costatenamet.value implicitnonlinearcontroltp1.value];
        data.description='dependent variables of the canonical system';

    case 'explicitnonlinearcontrol'
        % explicit nonlinear control variables for a specific arc
        if isfield(ocStruct.variable,'control')
            property=ocStruct.variable.control.(arcidentifier2field(arcidentifier)).property;
            data.value=ocStruct.variable.control.name([property.explicit{:}]==1 & [property.linear{:}]==0);
        else
            data.value=00;
        end
        data.type='cellchar';
        data.description='explicit nonlinear control variable name';

    case 'implicitnonlinearcontrol'
        % implicit nonlinear control variables for a specific arc
        property=ocStruct.variable.control.(arcidentifier2field(arcidentifier)).property;
        data.value='';
        data.type='cellchar';
        data.description='implicit nonlinear control variable name';

    case 'implicitnonlinearcontrolindex'
        % implicit nonlinear control variables for a specific arc
        try
            property=ocStruct.variable.control.(arcidentifier2field(arcidentifier)).property;
            data.value=find([property.explicit{:}]==0 & [property.linear{:}]==0);
        catch
            data.value=0;
        end
        data.type='integer';
        data.description='implicit nonlinear control variable indices';

    case 'linearcontrol'
        % linear control variables for a specific arc
        property=ocStruct.variable.control.(arcidentifier2field(arcidentifier)).property;
        data.value=ocStruct.variable.control.name([property.linear{:}]==1);
        data.type='cellchar';
        data.description='linear control variables';

    case 'maximizingvariable'
        % dependent variable names of the Hamiltonian maximizing condition
        controlname=retrieveimpulsemodelinformation(ocStruct,'controlname');
        lagrangemultipliercontrolname=retrieveimpulsemodelinformation(ocStruct,'lagrangemultipliercontrolname');
        data.value=[controlname.value lagrangemultipliercontrolname.value];
        data.type='cellchar';
        data.description='variable names of the Hamiltonian maximizing condition';

    case 'maximizingimplicitderivativevariable'
        % variable names which are used to derive the first order necessary
        % conditions of the Hamiltonian maximizing condition and are only
        % calculated implicitly
        try
            data=retrieveimpulsemodelinformation(ocStruct,'maximizingvariable');
            maximizingderivativevariable=retrieveimpulsemodelinformation(ocStruct,'maximizingderivativevariable',arcidentifier);
            zerolmmc=retrieveimpulsemodelinformation(ocStruct,'zerolmmc',arcidentifier);
            totalexplicitvariable=[maximizingderivativevariable.value zerolmmc.value];
            for ii=1:numel(totalexplicitvariable)
                data.value(strcmp(data.value,totalexplicitvariable{ii}))=[];
            end
        catch
            data.value=0;
        end
        data.description='implicit variable names for the first order necessary conditions of the Hamiltonian maximizing condition';
    case 'maximizingderivativevariable'
        % variable names which are used to derive the first order necessary
        % conditions of the Hamiltonian maximizing condition
        arc=arcidentifier2field(arcidentifier);
        if isfield(ocStruct,'foc') && isfield(ocStruct.foc,'generalinformation') && ...
                isfield(ocStruct.foc.generalinformation,arc) && ...
                isfield(ocStruct.foc.generalinformation.(arc),'maximizingderivativevariable')

            data.value=ocStruct.foc.generalinformation.(arc).maximizingderivativevariable.name;
        else
            explicitnonlinearcontrol=retrieveimpulsemodelinformation(ocStruct,'explicitnonlinearcontrol',arcidentifier);
            nonzerolmmc=retrieveimpulsemodelinformation(ocStruct,'nonzerolmmc',arcidentifier);
            data.value=[explicitnonlinearcontrol.value nonzerolmmc.value];
        end
        data.type='cellchar';
        data.description='variable names for the first order necessary conditions of the Hamiltonian maximizing condition';

    case 'maximizingexplicitvariable'
        % variable names which are used to derive the first order necessary
        % conditions of the Hamiltonian maximizing condition and are
        % calculated explicitly
        arc=arcidentifier2field(arcidentifier);
        if isfield(ocStruct,'foc') &&isfield(ocStruct.foc,'generalinformation') && ...
                isfield(ocStruct.foc.generalinformation,arc) && ...
                isfield(ocStruct.foc.generalinformation.(arc),'maximizingexplicitvariable')
            data.value=ocStruct.foc.generalinformation.(arc).maximizingexplicitvariable.name;
        else
            explicitnonlinearcontrol=retrieveimpulsemodelinformation(ocStruct,'explicitnonlinearcontrol',arcidentifier);
            nonzerolmmc=retrieveimpulsemodelinformation(ocStruct,'nonzerolmmc',arcidentifier);
            data.value=[explicitnonlinearcontrol.value nonzerolmmc.value];
        end
        data.type='cellchar';
        data.description='explicit variable names for the first order necessary conditions of the Hamiltonian maximizing condition';

    case 'maximizingimplicitvariable'
        % variable names which are implicitly calculated by the Hamiltonian
        % maximizing condition
        data=retrieveimpulsemodelinformation(ocStruct,'maximizingvariable');
        maximizingexplicitvariable=retrieveimpulsemodelinformation(ocStruct,'maximizingexplicitvariable',arcidentifier);
        zerolmmc=retrieveimpulsemodelinformation(ocStruct,'zerolmmc',arcidentifier);
        totalexplicitvariable=[maximizingexplicitvariable.value zerolmmc.value];
        for ii=1:numel(totalexplicitvariable)
            data.value(strcmp(data.value,totalexplicitvariable{ii}))=[];
        end
        data.description='implicit variable names of the Hamiltonian maximizing condition';


    case 'totalimplicitvariable'
        controlname=retrieveimpulsemodelinformation(ocStruct,'controlname');
        totalimplicitnonlinearcontrolindex=retrieveimpulsemodelinformation(ocStruct,'totalimplicitvariableindex');
        data.value=controlname.value(totalimplicitnonlinearcontrolindex.value);
        data.type='integer';
        data.description='index of the Lagrange multipliers corresponding to an inactive constraint';

    case 'totalimplicitvariableindex'
        arcidentifier=retrieveimpulsemodelinformation(ocStruct,'arcidentifier');
        totalimplicitnonlinearcontrolindex=[];
        for ii=1:numel(arcidentifier.value)
            implicitnonlinearcontrolindex=retrieveimpulsemodelinformation(ocStruct,'implicitnonlinearcontrolindex',arcidentifier.value{ii});
            totalimplicitnonlinearcontrolindex=[totalimplicitnonlinearcontrolindex implicitnonlinearcontrolindex.value];
        end
        data.value=unique(totalimplicitnonlinearcontrolindex);
        data.type='integer';
        data.description='index of the Lagrange multipliers corresponding to an inactive constraint';

    case 'equationvariablenameimplicit'
        % dependent variables of the canonical system for a specific arc
        data.type='cellchar';
        statename=retrieveimpulsemodelinformation(ocStruct,'statename');
        costatename=retrieveimpulsemodelinformation(ocStruct,'costatename');
        implicitcontrolname=retrieveimpulsemodelinformation(ocStruct,'maximizingimplicitvariable',arcidentifier);
        data.value=[statename.value costatename.value implicitcontrolname.value];
        data.description='dependent variables of the canonical system';

    case 'zerolmmcindex'
        % returns the index of the Lagrange multipliers corresponding to an
        % inactive constraint for a specific arc
        inequalitycontrolconstraintidentifier=retrieveimpulsemodelinformation(ocStruct,'inequalitycontrolconstraintidentifier');
        controlconstraintnum=retrieveimpulsemodelinformation(ocStruct,'inequalitycontrolconstraintnum');
        constraintcombination=retrieveimpulsemodelinformation(ocStruct,'constraintcombination',arcidentifier);
        zerolmmcindex=[];
        if ~strcmp(constraintcombination.value,'[]')
            for jj=1:numel(inequalitycontrolconstraintidentifier.value)
                testexpr=regexp(constraintcombination.value,['\<' inequalitycontrolconstraintidentifier.value{jj} '\>'],'ONCE');
                if isempty([testexpr{:}])
                    zerolmmcindex=[zerolmmcindex jj];
                end
            end
        else
            zerolmmcindex=1:controlconstraintnum.value;
        end
        data.value=zerolmmcindex;
        data.type='integer';
        data.description='index of the Lagrange multipliers corresponding to an inactive constraint';

    case 'zerolmmc'
        % returns the index of the Lagrange multipliers corresponding to an
        % inactive constraint for a specific arc
        lagrangemultipliercontrolname=retrieveimpulsemodelinformation(ocStruct,'lagrangemultipliercontrolname');
        zerolmmcindex=retrieveimpulsemodelinformation(ocStruct,'zerolmmcindex',arcidentifier);
        data.value=lagrangemultipliercontrolname.value(zerolmmcindex.value);
        data.type='cellchar';
        data.description='name of the Lagrange multipliers corresponding to an inactive constraint';

    case 'nonzerolmmc'
        % returns the name of the Lagrange multipliers corresponding to an
        % active constraint for a specific arc
        lagrangemultipliercontrolname=retrieveimpulsemodelinformation(ocStruct,'lagrangemultipliercontrolname');
        nonzerolmmcindex=retrieveimpulsemodelinformation(ocStruct,'nonzerolmmcindex',arcidentifier);
        data.value=lagrangemultipliercontrolname.value(nonzerolmmcindex.value);
        data.type='cellchar';
        data.description='name of the Lagrange multipliers corresponding to an active constraint';


    case 'nonzerolmmcindex'
        % returns the index of the Lagrange multipliers corresponding to an
        % active constraint for a specific arc
        inequalitycontrolconstraintnum=retrieveimpulsemodelinformation(ocStruct,'inequalitycontrolconstraintnum');
        zerolmmcindex=retrieveimpulsemodelinformation(ocStruct,'zerolmmcindex',arcidentifier);
        data.value=setdiff(1:inequalitycontrolconstraintnum.value,zerolmmcindex.value);
        data.type='integer';
        data.description='index of the Lagrange multipliers corresponding to an active constraint';

    case 'hamiltonianfunction'
        % general Hamilton function including the "Lagrange terms" from
        % inequality contraints
        if isfield(ocStruct,'hamiltonianfunction') && isfield(ocStruct.hamiltonianfunction,'term')
            data.value=ocStruct.hamiltonianfunction.term;
        else
            statenum=retrieveimpulsemodelinformation(ocStruct,'statenum');
            costatename=retrieveimpulsemodelinformation(ocStruct,'costatename');
            % integrand of the objective value function without discount factor
            g=retrieveimpulsemodelinformation(ocStruct,'objectiveintegrand');
            % state dynamics
            dxdt=retrieveimpulsemodelinformation(ocStruct,'statedynamics');
            data.value=g.value;
            for ii=1:statenum.value
                data.value=[data.value '+' costatename.value{ii} '*(' dxdt.value{ii} ')'];
            end
        end
        data.type='mathchar';
        data.description='Hamilton function P=g+lamba*f';

    case 'impulsehamiltonianfunction'
        % general Hamilton function including the "Lagrange terms" from
        % inequality contraints
        if isfield(ocStruct,'impulsehamiltonianfunction') && isfield(ocStruct.impulsehamiltonianfunction,'term')
            data.value=ocStruct.hamiltonianfunction.term;
        else
            statenum=retrieveimpulsemodelinformation(ocStruct,'statenum');
            costatename=retrieveimpulsemodelinformation(ocStruct,'costatenamer');
            % integrand of the objective value function without discount factor
            g=retrieveimpulsemodelinformation(ocStruct,'objectivesummand');
            % state dynamics
            stateevent=retrieveimpulsemodelinformation(ocStruct,'stateevent');
            data.value=g.value;
            for ii=1:statenum.value
                data.value=[data.value '+' costatename.value{ii} '*(' stateevent.value{ii} ')'];
            end
        end
        data.type='mathchar';
        data.description='Hamilton function P=g+lamba*f';


    case 'pontryaginfunction'
        % general Hamilton function including the "Lagrange terms" from
        % inequality contraints
        if isfield(ocStruct,'pontryaginfunction') && isfield(ocStruct.pontryaginfunction,'term')
            data.value=ocStruct.pontryaginfunction.term;
        else
            statenum=retrieveimpulsemodelinformation(ocStruct,'statenum');
            costatename=retrieveimpulsemodelinformation(ocStruct,'costatename');
            inequalitycontrolconstraint=retrieveimpulsemodelinformation(ocStruct,'inequalitycontrolconstraint');
            inequalitycontrolconstraintnum=retrieveimpulsemodelinformation(ocStruct,'inequalitycontrolconstraintnum');
            lmmc=retrieveimpulsemodelinformation(ocStruct,'lagrangemultipliercontrolname');
            % integrand of the objective value function without discount factor
            g=retrieveimpulsemodelinformation(ocStruct,'objectiveintegrand');
            % state dynamics
            dxdt=retrieveimpulsemodelinformation(ocStruct,'statedynamics');
            data.value=g.value;
            for ii=1:statenum.value
                data.value=[data.value '+' costatename.value{ii} '*(' dxdt.value{ii} ')'];
            end
            for ii=1:inequalitycontrolconstraintnum.value
                data.value=[data.value '+' lmmc.value{ii} '*(' inequalitycontrolconstraint.value{ii} ')'];
            end
        end
        data.type='mathchar';
        data.description='general Hamilton function P=H+mu*mc+nu*sc';

    case 'pontryaginfunctionDx'
        % derivative of the Pontryaginfunction with respect to the state
        % variable(s)
        if isfield(ocStruct,'pontryaginfunction') && isfield(ocStruct.pontryaginfunction,'derivative') ...
                && isfield(ocStruct.pontryaginfunction.derivative,'Dx')
            data.value=ocStruct.pontryaginfunction.derivative.Dx.term;
        else
            statename=retrieveimpulsemodelinformation(ocStruct,'statename');
            statevector=cell2vectorstring(statename.value);
            pontryaginfunction=retrieveimpulsemodelinformation(ocStruct,'pontryaginfunction');
            data.value=string2cell(ocmatjacobian(['[' pontryaginfunction.value ']'],statevector,symkernel),'vector');
        end
        data.type='mathcellchar';
        data.description='derivative of the Pontryaginfunction with respect to the state';

    case 'pontryaginfunctionDX'
        % derivative of the Pontryaginfunction with respect to the state
        % AND costate variable(s)
        costatename=retrieveimpulsemodelinformation(ocStruct,'costatename');
        statename=retrieveimpulsemodelinformation(ocStruct,'statename');
        statevector=cell2vectorstring([statename.value costatename.value]);
        pontryaginfunction=retrieveimpulsemodelinformation(ocStruct,'specificpontryaginfunction',arcidentifier,symkernel);
        data.value=string2cell(ocmatjacobian(['[' pontryaginfunction.value ']'],statevector,symkernel),'vector');
        data.type='mathcellchar';
        data.description='derivative of the Pontryaginfunction with respect to the state and costate';


    case 'pontryaginfunctionDu'
        % derivative of the Pontryaginfunction with respect to the control
        % variable(s)
        if isfield(ocStruct,'pontryaginfunction') && isfield(ocStruct.pontryaginfunction,'derivative') ...
                && isfield(ocStruct.pontryaginfunction.derivative,'Du')
            data.value=ocStruct.pontryaginfunction.derivative.Du.term;
        else
            controlname=retrieveimpulsemodelinformation(ocStruct,'controlname');
            if ~isempty(controlname.value)
                controlvector=cell2vectorstring(controlname.value);
                pontryaginfunction=retrieveimpulsemodelinformation(ocStruct,'pontryaginfunction');
                data.value=string2cell(ocmatjacobian(['[' pontryaginfunction.value ']'],controlvector,symkernel),'vector');
            else
                data.value=[];
            end
        end
        data.type='mathcellchar';
        data.description='derivative of the Pontryaginfunction with respect to the control';

    case 'pontryaginfunctionDx2'
        % second order derivative of the Pontryaginfunction with respect to
        % the state variable(s)
        if isfield(ocStruct,'pontryaginfunction') && isfield(ocStruct.pontryaginfunction,'derivative') ...
                && isfield(ocStruct.pontryaginfunction.derivative,'Dx2')
            data.value=ocStruct.pontryaginfunction.derivative.Dx2.term;
        else
            statename=retrieveimpulsemodelinformation(ocStruct,'statename');
            statevector=cell2vectorstring(statename.value);
            pontryaginfunction=retrieveimpulsemodelinformation(ocStruct,'pontryaginfunction');
            data.value=string2cell(ocmathessian(pontryaginfunction.value,statevector,symkernel),'matrix');
        end
        data.type='mathcellchar';
        data.description='second order derivative of the Pontryaginfunction with respect to the state';

    case 'pontryaginfunctionDu2'
        % second order derivative of the Pontryaginfunction with respect to
        % the control variable(s)
        if isfield(ocStruct,'pontryaginfunction') && isfield(ocStruct.pontryaginfunction,'derivative') ...
                && isfield(ocStruct.pontryaginfunction.derivative,'Du2')
            data.value=ocStruct.pontryaginfunction.derivative.Du2.term;
        else
            controlname=retrieveimpulsemodelinformation(ocStruct,'controlname');
            controlvector=cell2vectorstring(controlname.value);
            pontryaginfunction=retrieveimpulsemodelinformation(ocStruct,'pontryaginfunction');
            if ~isempty(controlname.value)
                data.value=string2cell(ocmathessian(pontryaginfunction.value,controlvector,symkernel),'matrix');
            else
                data.value=[];
            end
        end
        data.type='mathcellchar';
        data.description='second order derivative of the Pontryaginfunction with respect to the control';

    case {'pontryaginfunctionDuDx','pontryaginfunctionDxDu'}
        % second order derivative of the Pontryaginfunction with respect to
        % the control variable(s)
        if isfield(ocStruct,'pontryaginfunction') && isfield(ocStruct.pontryaginfunction,'derivative') ...
                && isfield(ocStruct.pontryaginfunction.derivative,'DuDx')
            data.value=ocStruct.pontryaginfunction.derivative.DuDx.term;
        else
            controlname=retrieveimpulsemodelinformation(ocStruct,'controlname');
            if ~isempty(controlname.value)
                controlvector=cell2vectorstring(controlname.value);
                statename=retrieveimpulsemodelinformation(ocStruct,'statename');
                statevector=cell2vectorstring(statename.value);
                pontryaginfunction=retrieveimpulsemodelinformation(ocStruct,'pontryaginfunction');
                J=removematrixstring(ocmatjacobian(['[' pontryaginfunction.value ']'],controlvector,symkernel));
                data.value=string2cell(ocmatjacobian(J,statevector,symkernel),'matrix');
            else
                data.value=[];
            end
        end
        data.type='mathcellchar';
        data.description='second order derivative of the Pontryaginfunction with respect to the control';

    case 'impulsepontryaginfunction'
        % general Hamilton function including the "Lagrange terms" from
        % inequality contraints
        if isfield(ocStruct,'impulsepontryaginfunction') && isfield(ocStruct.impulsepontryaginfunction,'term')
            data.value=ocStruct.impulsepontryaginfunction.term;
        else
            statenum=retrieveimpulsemodelinformation(ocStruct,'statenum');
            costatename=retrieveimpulsemodelinformation(ocStruct,'costatenamer');
            % integrand of the objective value function without discount factor
            g=retrieveimpulsemodelinformation(ocStruct,'objectivesummand');
            % state dynamics
            stateevent=retrieveimpulsemodelinformation(ocStruct,'stateevent');
            data.value=g.value;
            for ii=1:statenum.value
                data.value=[data.value '+' costatename.value{ii} '*(' stateevent.value{ii} ')'];
            end
        end
        data.type='mathchar';
        data.description='general Impulse Hamilton function P=H+mu*mc+nu*sc';


    case 'impulsepontryaginfunctionDx'
        % derivative of the Pontryaginfunction with respect to the state
        % variable(s)
        if isfield(ocStruct,'impulsepontryaginfunction') && isfield(ocStruct.impulsepontryaginfunction,'derivative') ...
                && isfield(ocStruct.impulsepontryaginfunction.derivative,'Dx')
            data.value=ocStruct.impulsepontryaginfunction.derivative.Dx.term;
        else
            statename=retrieveimpulsemodelinformation(ocStruct,'statenamel');
            statevector=cell2vectorstring(statename.value);
            impulsepontryaginfunction=retrieveimpulsemodelinformation(ocStruct,'impulsepontryaginfunction');
            data.value=string2cell(ocmatjacobian(['[' impulsepontryaginfunction.value ']'],statevector,symkernel),'vector');
        end
        data.type='mathcellchar';
        data.description='derivative of the impulse Pontryaginfunction with respect to the state';


    case 'impulsepontryaginfunctionDu'
        % derivative of the Pontryaginfunction with respect to the control
        % variable(s)
        if isfield(ocStruct,'impulsepontryaginfunction') && isfield(ocStruct.impulsepontryaginfunction,'derivative') ...
                && isfield(ocStruct.impulsepontryaginfunction.derivative,'Du')
            data.value=ocStruct.impulsepontryaginfunction.derivative.Du.term;
        else
            icontrolname=retrieveimpulsemodelinformation(ocStruct,'icontrolname');
            if ~isempty(icontrolname.value)
                icontrolvector=cell2vectorstring(icontrolname.value);
                impulsepontryaginfunction=retrieveimpulsemodelinformation(ocStruct,'impulsepontryaginfunction');
                data.value=string2cell(ocmatjacobian(['[' impulsepontryaginfunction.value ']'],icontrolvector,symkernel),'vector');
            else
                data.value=[];
            end
        end
        data.type='mathcellchar';
        data.description='derivative of the impulse Pontryaginfunction with respect to the control';

    case 'impulsepontryaginfunctionDtau'
        % derivative of the Pontryaginfunction with respect to the control
        % variable(s)
        if isfield(ocStruct,'impulsepontryaginfunction') && isfield(ocStruct.impulsepontryaginfunction,'derivative') ...
                && isfield(ocStruct.impulsepontryaginfunction.derivative,'Dtau')
            data.value=ocStruct.impulsepontryaginfunction.derivative.Dtau.term;
        else
            impulsetime=retrieveimpulsemodelinformation(ocStruct,'impulsetime');
            if ~isempty(impulsetime.value)
                impulsetimevector=cell2vectorstring(impulsetime.value);
                impulsepontryaginfunction=retrieveimpulsemodelinformation(ocStruct,'impulsepontryaginfunction');
                data.value=string2cell(ocmatjacobian(['[' impulsepontryaginfunction.value ']'],impulsetimevector,symkernel),'vector');
            else
                data.value=[];
            end
        end
        data.type='mathcellchar';
        data.description='derivative of the impulse Pontryaginfunction with respect to the impulse time';

    case 'DsalvagevalueDx'
        % derivative of the Pontryaginfunction with respect to the state
        % variable(s)
        statename=retrieveimpulsemodelinformation(ocStruct,'statenamer');
        statevector=cell2vectorstring(statename.value);
        salvagevalue=retrieveimpulsemodelinformation(ocStruct,'salvagevalue');
        data.value=string2cell(ocmatjacobian(['[' salvagevalue.value ']'],statevector,symkernel),'vector');
        data.type='mathcellchar';
        data.description='derivative of the Pontryaginfunction with respect to the state';

    case 'DsalvagevalueDT'
        % derivative of the Pontryaginfunction with respect to the state
        % variable(s)
        endtime=retrieveimpulsemodelinformation(ocStruct,'endtime');
        salvagevalue=retrieveimpulsemodelinformation(ocStruct,'salvagevalue');
        endtimevector=cell2vectorstring(endtime.value);
        data.value=string2cell(ocmatjacobian(['[' salvagevalue.value ']'],endtimevector,symkernel),'vector');
        data.type='mathcellchar';
        data.description='derivative of the Pontryaginfunction with respect to the state';

    case 'optimalcontroldynamicsleftside'
        controlname=retrieveimpulsemodelinformation(ocStruct,'controlname');
        controlvector=cell2vectorstring(controlname.value);
        zerolmmcindex=retrieveimpulsemodelinformation(ocStruct,'zerolmmcindex',arcidentifier);
        lagrangemultipliercontrolname=retrieveimpulsemodelinformation(ocStruct,'lagrangemultipliercontrolname');
        lagrangemultipliercontrolname.value(zerolmmcindex.value)=[];
        lagrangemultipliervector=cell2vectorstring(lagrangemultipliercontrolname.value);
        pontryaginfunction=retrieveimpulsemodelinformation(ocStruct,'pontryaginfunction');
        pontryaginfunction4arc=retrieveimpulsemodelinformation(ocStruct,'pontryaginfunction4arc',arcidentifier,symkernel);
        pontryaginfunctionDu2=sym(ocmathessian(pontryaginfunction.value,controlvector,symkernel));
        if ~isempty(lagrangemultipliervector)
            pontryaginfunctionDlmmcDu=sym(ocmatjacobian(removematrixstring(ocmatjacobian(['[' pontryaginfunction4arc.value ']'],controlvector,symkernel)),lagrangemultipliervector,symkernel));
        else
            pontryaginfunctionDlmmcDu=sym([]);
        end
        if ~isempty(lagrangemultipliercontrolname.value)
            data.value=string2cell(char([pontryaginfunctionDu2 pontryaginfunctionDlmmcDu;pontryaginfunctionDlmmcDu.' sym(zeros(length(lagrangemultipliercontrolname.value)))]),'charmatrix');
        else
            if numel(pontryaginfunctionDu2)>1
                data.value=string2cell(char(pontryaginfunctionDu2),'charmatrix');
            else
                data.value{1}=char(pontryaginfunctionDu2);
            end
        end
        data.type='mathcellchar';
        data.description='left side matrix of the optimal control dynamics';


    case {'pontryaginfunctionDuDX','pontryaginfunctionDXDu'}
        try
            data.value=ocStruct.pontryaginfunction.derivative.DuDX.term;
        catch
            controlname=retrieveimpulsemodelinformation(ocStruct,'controlname');
            controlvector=cell2vectorstring(controlname.value);
            statename=retrieveimpulsemodelinformation(ocStruct,'statename');
            costatename=retrieveimpulsemodelinformation(ocStruct,'costatename');
            statevector=cell2vectorstring([statename.value costatename.value]);
            pontryaginfunction=retrieveimpulsemodelinformation(ocStruct,'pontryaginfunction');
            J=removematrixstring(ocmatjacobian(['[' pontryaginfunction.value ']'],controlvector,symkernel));
            data.value=string2cell(ocmatjacobian(J,statevector,symkernel),'matrix');
        end
        data.type='mathcellchar';
        data.description='right side of the optimal control dynamics (-Pux*dotx-Pul*dotl)';

    case {'pontryaginfunctionDuDlmmc','pontryaginfunctionDlmmcDu'}
        try
            data.value=ocStruct.pontryaginfunction.derivative.DuDlmmc.term;
        catch
            controlname=retrieveimpulsemodelinformation(ocStruct,'controlname');
            controlvector=cell2vectorstring(controlname.value);
            lagrangemultipliercontrolname=retrieveimpulsemodelinformation(ocStruct,'lagrangemultipliercontrolname');
            lagrangemultipliervector=cell2vectorstring(lagrangemultipliercontrolname.value);
            pontryaginfunction=retrieveimpulsemodelinformation(ocStruct,'pontryaginfunction');
            J=removematrixstring(ocmatjacobian(['[' pontryaginfunction.value ']'],controlvector,symkernel));
            if ~isempty(lagrangemultipliervector)
                data.value=string2cell(ocmatjacobian(J,lagrangemultipliervector,symkernel),'matrix');
            else
                data.value='';
            end
        end
        data.type='mathcellchar';
        data.description='right side of the optimal control dynamics (-Pux*dotx-Pul*dotl)';

    case {'pontryaginfunctionDXDlmmc','pontryaginfunctionDlmmcDX'}
        statename=retrieveimpulsemodelinformation(ocStruct,'statename');
        costatename=retrieveimpulsemodelinformation(ocStruct,'costatename');
        statevector=cell2vectorstring([statename.value costatename.value]);
        zerolmmcindex=retrieveimpulsemodelinformation(ocStruct,'zerolmmcindex',arcidentifier);
        lagrangemultipliercontrolname=retrieveimpulsemodelinformation(ocStruct,'lagrangemultipliercontrolname');
        lagrangemultipliercontrolname.value(zerolmmcindex.value)=[];
        lagrangemultipliervector=cell2vectorstring(lagrangemultipliercontrolname.value);
        pontryaginfunction=retrieveimpulsemodelinformation(ocStruct,'pontryaginfunction');
        J=removematrixstring(ocmatjacobian(['[' pontryaginfunction.value ']'],statevector,symkernel));
        if ~isempty(lagrangemultipliervector)
            data.value=string2cell(ocmatjacobian(J,lagrangemultipliervector,symkernel),'matrix');
        else
            data.value='';
        end
        data.type='mathcellchar';
        data.description='right side of the optimal control dynamics (-Pux*dotx-Pul*dotl)';

    case {'specificpontryaginfunctionDuDlmmc','specificpontryaginfunctionDlmmcDu'}
        try
            pontryaginfunctionDuDlmmc=mycell2sym(ocStruct.pontryaginfunction.derivative.DuDlmmc.term,'matrix');
        catch
            controlname=retrieveimpulsemodelinformation(ocStruct,'controlname');
            controlvector=cell2vectorstring(controlname.value);
            lagrangemultipliercontrolname=retrieveimpulsemodelinformation(ocStruct,'lagrangemultipliercontrolname');
            lagrangemultipliervector=cell2vectorstring(lagrangemultipliercontrolname.value);
            pontryaginfunction=retrieveimpulsemodelinformation(ocStruct,'pontryaginfunction');
            pontryaginfunctionDuDlmmc=removematrixstring(ocmatjacobian(['[' pontryaginfunction.value ']'],controlvector,symkernel));
            pontryaginfunctionDuDlmmc=ocmatjacobian(pontryaginfunctionDuDlmmc,lagrangemultipliervector,symkernel);
        end
        zerolmmcindex=retrieveimpulsemodelinformation(ocStruct,'zerolmmcindex',arcidentifier);
        pontryaginfunctionDuDlmmc(:,zerolmmcindex.value)=[];
        data.value=string2cell(removematrixstring(char(pontryaginfunctionDuDlmmc)),'charmatrix');
        data.type='mathcellchar';
        data.description='right side of the optimal control dynamics (-Pux*dotx-Pul*dotl)';

    case {'zeros4optimalcontroldynamics'}
        zerolmmcindex=retrieveimpulsemodelinformation(ocStruct,'zerolmmcindex',arcidentifier);
        lagrangemultipliercontrolname=retrieveimpulsemodelinformation(ocStruct,'lagrangemultipliercontrolname');
        lagrangemultipliercontrolname.value(zerolmmcindex.value)=[];
        l=length(lagrangemultipliercontrolname.value);
        if l>0
            data.value=string2cell(removematrixstring(char(sym(zeros(l)))),'charmatrix');
        else
            data.value='';
        end
        data.type='mathcellchar';
        data.description='right side of the optimal control dynamics (-Pux*dotx-Pul*dotl)';

    case 'optimalcontroldynamicsindex'
        try
            data.value=ocStruct.foc.adjointsystem.optimalcontroldynamics.(arcidentifier2field(arcidentifier)).leftsideindex.value;
        catch
            maximizingimplicitvariable=retrieveimpulsemodelinformation(ocStruct,'maximizingimplicitvariable',arcidentifier);
            controlname=retrieveimpulsemodelinformation(ocStruct,'controlname');
            zerolmmcindex=retrieveimpulsemodelinformation(ocStruct,'zerolmmcindex',arcidentifier);
            lagrangemultipliercontrolname=retrieveimpulsemodelinformation(ocStruct,'lagrangemultipliercontrolname');
            lagrangemultipliercontrolname.value(zerolmmcindex.value)=[];
            controldynamicsstates=[controlname.value lagrangemultipliercontrolname.value];
            counter=0;
            data.value=[];
            for ii=1:length(maximizingimplicitvariable.value)
                idx=strfind(controldynamicsstates,maximizingimplicitvariable.value{ii});
                if ~isempty([idx{:}])
                    counter=counter+1;
                    for jj=1:length(controldynamicsstates)
                        if ~isempty(idx{jj})
                            data.value(counter)=jj;
                            break
                        end
                    end

                end
            end
        end
        data.type='vector';
        data.description='left side matrix of the optimal control dynamics';

    case 'pontryaginfunction4arc'
        % returns the specific Pontryagin functions for every considered
        % constraint combination. It is assumed that the constraints are
        % given by ODEs (statedynamics) and inequality constraints
        % (control-state constraints)
        pontryaginfunction=retrieveimpulsemodelinformation(ocStruct,'pontryaginfunction');
        zerolmmcindex=retrieveimpulsemodelinformation(ocStruct,'zerolmmcindex',arcidentifier);
        lagrangemultipliercontrolname=retrieveimpulsemodelinformation(ocStruct,'lagrangemultipliercontrolname');
        lmmcequalzero='';
        ii=0;
        for idx=zerolmmcindex.value
            ii=ii+1;
            lmmcequalzero{ii}=[lagrangemultipliercontrolname.value{idx} '=0'];
        end
        zerolmmcvector=cell2vectorstring(lmmcequalzero);
        if ~isempty(zerolmmcvector)
            pontryaginfunction.value=ocmatsubs(pontryaginfunction.value,zerolmmcvector,symkernel);
        end
        data.value=pontryaginfunction.value;
        data.type='mathchar';
        data.description='specific Hamilton function for actual constraint combination';

    case 'impulsepontryaginfunction4arc'
        % returns the specific Pontryagin functions for every considered
        % constraint combination. It is assumed that the constraints are
        % given by ODEs (statedynamics) and inequality constraints
        % (control-state constraints)
        impulsepontryaginfunction=retrieveimpulsemodelinformation(ocStruct,'impulsepontryaginfunction');
        zerolmmcindex=retrieveimpulsemodelinformation(ocStruct,'zerolmmcindex',arcidentifier);
        lagrangemultipliercontrolname=retrieveimpulsemodelinformation(ocStruct,'lagrangemultipliercontrolname');
        lmmcequalzero='';
        ii=0;
        for idx=zerolmmcindex.value
            ii=ii+1;
            lmmcequalzero{ii}=[lagrangemultipliercontrolname.value{idx} '=0'];
        end
        zerolmmcvector=cell2vectorstring(lmmcequalzero);
        if ~isempty(zerolmmcvector)
            impulsepontryaginfunction.value=ocmatsubs(impulsepontryaginfunction.value,zerolmmcvector,symkernel);
        end
        data.value=impulsepontryaginfunction.value;
        data.type='mathchar';
        data.description='specific Hamilton function for actual constraint combination';

    case 'impulsemaximizingsolution'
        % Pontryagin's maximumprinciple is used to derive explicit formulas
        % of control and Lagrange multiplier values satisfying the first
        % order necessary conditions
        solution=[];
        icontrolname=retrieveimpulsemodelinformation(ocStruct,'icontrolname');
        icontrolnum=retrieveimpulsemodelinformation(ocStruct,'icontrolnum');
        impulsepontryaginfunction4arc=retrieveimpulsemodelinformation(ocStruct,'impulsepontryaginfunction4arc',arcidentifier,symkernel);
        equationvariable=cell2vectorstring([icontrolname.value]);
        derivativevariable=cell2vectorstring(icontrolname.value);
        equation=removematrixstring(ocmatjacobian(['[' impulsepontryaginfunction4arc.value ']'],derivativevariable,symkernel));
        % ocmatsolve is an adaptation of the native MATLAB command
        % solve for the symbolic toolbox relying
        if ~strcmp(equation,'[]')
            solution=ocmatsolve(equation,equationvariable,symkernel);
            if isempty(solution)
                %ocmatmsg('Explicit solution could not be found.\n')
                return
            end
        end
        solution=orderfields(solution,icontrolname.value);
        for ii=1:numel(solution)
            for jj=1:icontrolnum.value
                solution(ii).icontrol{jj}=solution(ii).(icontrolname.value{jj});
            end
        end
        data.value=solution;
        data.type='struct';
        data.description='solution(s) derived from the foc of the Hamiltonian maximizing condition';

    case 'maximizingsolution'
        % Pontryagin's maximumprinciple is used to derive explicit formulas
        % of control and Lagrange multiplier values satisfying the first
        % order necessary conditions
        solution=[];
        totalmaximizingvariable=retrieveimpulsemodelinformation(ocStruct,'maximizingvariable');
        if ~isempty(totalmaximizingvariable.value)
            maximizingexplicitvariable=retrieveimpulsemodelinformation(ocStruct,'maximizingexplicitvariable',arcidentifier);
            maximizingimplicitvariable=retrieveimpulsemodelinformation(ocStruct,'maximizingimplicitvariable',arcidentifier);
            zerolmmc=retrieveimpulsemodelinformation(ocStruct,'zerolmmc',arcidentifier);
            nonzerolmmc=retrieveimpulsemodelinformation(ocStruct,'nonzerolmmc',arcidentifier);
            controlname=retrieveimpulsemodelinformation(ocStruct,'controlname');
            lagrangemultipliercontrolname=retrieveimpulsemodelinformation(ocStruct,'lagrangemultipliercontrolname');
            controlnum=retrieveimpulsemodelinformation(ocStruct,'controlnum');
            inequalitycontrolconstraintnum=retrieveimpulsemodelinformation(ocStruct,'inequalitycontrolconstraintnum');
            pontryaginfunction4arc=retrieveimpulsemodelinformation(ocStruct,'pontryaginfunction4arc',arcidentifier,symkernel);
            idx1=ocmatrank(pontryaginfunction4arc.value,cell2vectorstring(nonzerolmmc.value),cell2vectorstring(maximizingimplicitvariable.value));
            if ~isempty(idx1)
                nonzerolmmc.value(idx1)=[];
                maximizingvariable.value=[controlname.value nonzerolmmc.value];
            else
                maximizingvariable.value=maximizingexplicitvariable.value;
            end
            idx=ocmatrank(pontryaginfunction4arc.value,cell2vectorstring(maximizingvariable.value),cell2vectorstring(maximizingexplicitvariable.value));

            if ~isempty(maximizingexplicitvariable.value)
                equationvariable=cell2vectorstring([maximizingexplicitvariable.value]);
                derivativevariable=cell2vectorstring([maximizingvariable.value(idx)]);
                equation=removematrixstring(ocmatjacobian(['[' pontryaginfunction4arc.value ']'],derivativevariable,symkernel));
                % ocmatsolve is an adaptation of the native MATLAB command
                % solve for the symbolic toolbox relying
                if ~strcmp(equation,'[]')
                    solution=ocmatsolve(equation,equationvariable,symkernel);
                    if isempty(solution)
                        %ocmatmsg('Explicit solution could not be found.\n')
                        return
                    end
                end
                if ~isempty(solution)
                    for jj=1:numel(solution)
                        % the implicit control values are added to the solution
                        % structure, where the value is the variable itself
                        for ii=1:numel(maximizingimplicitvariable.value)
                            solution(jj).(maximizingimplicitvariable.value{ii})=maximizingimplicitvariable.value{ii};
                        end
                    end
                    for jj=1:numel(solution)
                        for ii=1:numel(zerolmmc.value)
                            % the Lagrange multipliers being zero (inactive
                            % constraint) are added
                            solution(jj).(zerolmmc.value{ii})='0';
                        end
                    end
                else
                    % the implicit control values are added to the solution
                    % structure, where the value is the variable itself
                    for ii=1:numel(maximizingimplicitvariable.value)
                        solution.(maximizingimplicitvariable.value{ii})=maximizingimplicitvariable.value{ii};
                    end
                    for ii=1:numel(zerolmmc.value)
                        % the Lagrange multipliers being zero (inactive
                        % constraint) are added
                        solution.(zerolmmc.value{ii})='0';
                    end
                end
                % the order of the solution fields are controlname,
                % lmmcname
            elseif ~isempty(maximizingimplicitvariable.value)
                for ii=1:numel(maximizingimplicitvariable.value)
                    solution.(maximizingimplicitvariable.value{ii})=maximizingimplicitvariable.value{ii};
                end
                for ii=1:numel(zerolmmc.value)
                    % the Lagrange multipliers being zero (inactive
                    % constraint) are added
                    solution.(zerolmmc.value{ii})='0';
                end
            else
                solution=[];
            end
            solution=orderfields(solution,totalmaximizingvariable.value);
            for ii=1:numel(solution)
                for jj=1:controlnum.value
                    solution(ii).control{jj}=solution(ii).(controlname.value{jj});
                end
                for jj=1:inequalitycontrolconstraintnum.value
                    solution(ii).lagrangemultcc{jj}=solution(ii).(lagrangemultipliercontrolname.value{jj});
                end
            end
        else
            solution=[];
        end
        data.value=solution;
        data.type='struct';
        data.description='solution(s) derived from the foc of the Hamiltonian maximizing condition';



    case 'adjointsystem'
        % the adjoint system, where each cell contains one costate dynamics
        arc=arcidentifier2field(arcidentifier);
        if isfield(ocStruct,'foc') && isfield(ocStruct.foc,'adjointsystem') ...
                && isfield(ocStruct.foc.adjointsystem,'dynamics') ...
                && isfield(ocStruct.foc.adjointsystem.dynamics,arc) ...
                && isfield(ocStruct.foc.adjointsystem.dynamics.(arc),'ode')
            data.value=ocStruct.foc.adjointsystem.dynamics.(arc).ode.term;
        else
            pontryaginfunctionDx=retrieveimpulsemodelinformation(ocStruct,'pontryaginfunctionDx','',symkernel);
            statenum=retrieveimpulsemodelinformation(ocStruct,'statenum');
            costatename=retrieveimpulsemodelinformation(ocStruct,'costatename');
            discountratevariable=retrieveimpulsemodelinformation(ocStruct,'discountratevariable');
            for ii=1:statenum.value
                data.value{ii}=[discountratevariable.value '*' costatename.value{ii} '-(' pontryaginfunctionDx.value{ii} ')'];
            end
        end
        data.type='mathcellchar';
        data.description='costate dynamics';

    case 'adjointevent'
        % the adjoint system, where each cell contains one costate dynamics
        arc=arcidentifier2field(arcidentifier);
        if isfield(ocStruct,'foc') && isfield(ocStruct.foc,'adjointsystem') ...
                && isfield(ocStruct.foc.adjointsystem,'dynamics') ...
                && isfield(ocStruct.foc.adjointsystem.dynamics,arc) ...
                && isfield(ocStruct.foc.adjointsystem.dynamics.(arc),'evt')
            data.value=ocStruct.foc.adjointsystem.dynamics.(arc).evt.term;
        else
            impulsepontryaginfunctionDx=retrieveimpulsemodelinformation(ocStruct,'impulsepontryaginfunctionDx','',symkernel);
            statenum=retrieveimpulsemodelinformation(ocStruct,'statenum');
            for ii=1:statenum.value
                data.value{ii}=['-(' impulsepontryaginfunctionDx.value{ii} ')'];
            end
        end
        data.type='mathcellchar';
        data.description='costate event';

    case 'specificstatedynamics'
        % the canonical system consisting of the state dynamics, adjoint
        % system and possible algebraic equations for implicitly given
        % control values.
        statedynamics=retrieveimpulsemodelinformation(ocStruct,'statedynamics');
        optimalvalue=retrieveimpulsemodelinformation(ocStruct,'optimalvalue',arcidentifier);
        maximizingvariable=retrieveimpulsemodelinformation(ocStruct,'maximizingvariable');
        dxdt=statedynamics.value;

        if ~isempty(optimalvalue.value)
            for ii=1:numel(maximizingvariable.value)
                for jj=1:numel(dxdt)
                    dxdt{jj}=ocmatsubs(dxdt{jj},[maximizingvariable.value{ii} '=' optimalvalue.value.(maximizingvariable.value{ii})],symkernel);
                end
            end
        end
        data.value=dxdt;
        data.type='mathcellchar';
        data.description='state dynamics';

    case 'specificadjointsystem'
        % the canonical system consisting of the state dynamics, adjoint
        % system and possible algebraic equations for implicitly given
        % control values.
        adjointsystem=retrieveimpulsemodelinformation(ocStruct,'adjointsystem','',symkernel);
        optimalvalue=retrieveimpulsemodelinformation(ocStruct,'optimalvalue',arcidentifier);
        maximizingvariable=retrieveimpulsemodelinformation(ocStruct,'maximizingvariable');
        dldt=adjointsystem.value;

        if ~isempty(optimalvalue.value)
            for ii=1:numel(maximizingvariable.value)
                for jj=1:numel(dldt)
                    dldt{jj}=ocmatsubs(dldt{jj},[maximizingvariable.value{ii} '=' optimalvalue.value.(maximizingvariable.value{ii})]);
                end
            end
        end
        data.value=dldt;
        data.type='mathcellchar';
        data.description='costate dynamics';

    case 'adjointstateexplicit'
        % the canonical system consisting of the state dynamics, adjoint
        % system and possible algebraic equations for implicitly given
        % control values.
        adjointsystem=retrieveimpulsemodelinformation(ocStruct,'adjointsystem','',symkernel);
        optimalvalue=retrieveimpulsemodelinformation(ocStruct,'optimalvalue',arcidentifier);
        maximizingvariable=retrieveimpulsemodelinformation(ocStruct,'maximizingvariable');
        costatenametp1=retrieveimpulsemodelinformation(ocStruct,'costatenametp1');
        statenum=retrieveimpulsemodelinformation(ocStruct,'statenum');
        for ii=1:numel(maximizingvariable.value)
            for jj=1:numel(adjointsystem.value)
                adjointsystem.value{jj}=ocmatsubs(adjointsystem.value{jj},[maximizingvariable.value{ii} '=' optimalvalue.value.(maximizingvariable.value{ii})]);
            end
        end
        equationvariable=cell2vectorstring(costatenametp1.value);
        equation=cell2vectorstring(adjointsystem.value);
        % ocmatsolve is an adaptation of the native MATLAB command
        % solve for the symbolic toolbox relying
        solution=ocmatsolve(equation,equationvariable,symkernel);
        if ~isempty(solution) && numel(solution)==1
            for ii=1:statenum.value
                data.value{ii}=ocmatsimple(solution.(costatenametp1.value{ii}),symkernel);
            end
        else
            data.value=[];
        end
        data.type='mathcellchar';
        data.description='state dynamics, adjoint system (and algebraic equations for implicit control values)';

    case 'specificcanonicalsystem'
        % the canonical system consisting of the state dynamics, adjoint
        % system and possible algebraic equations for implicitly given
        % control values.
        adjointsystem=retrieveimpulsemodelinformation(ocStruct,'adjointsystem','',symkernel);
        statedynamics=retrieveimpulsemodelinformation(ocStruct,'statedynamics');
        optimalvalue=retrieveimpulsemodelinformation(ocStruct,'optimalvalue',arcidentifier);
        maximizingvariable=retrieveimpulsemodelinformation(ocStruct,'maximizingvariable');
        dxdt=[statedynamics.value adjointsystem.value];
        for ii=1:numel(maximizingvariable.value)
            for jj=1:numel(dxdt)
                if ~(dxdt{jj}==sym('0'))
                    dxdt{jj}=ocmatsubs(dxdt{jj},[maximizingvariable.value{ii} '=' optimalvalue.value.(maximizingvariable.value{ii})],symkernel);
                end
            end
        end
        data.value=dxdt;
        data.type='mathcellchar';
        data.description='state dynamics, adjoint system (and algebraic equations for implicit control values)';

    case 'equilibriumequation'
        specificcanonicalsystem=retrieveimpulsemodelinformation(ocStruct,'specificcanonicalsystem',arcidentifier,symkernel);
        data.value=specificcanonicalsystem.value;
        data.type='mathcellchar';
        data.description='';

    case 'specificcanonicalsystemexplicit'
        % the canonical system consisting of the state dynamics, adjoint
        % system and possible algebraic equations for implicitly given
        % control values.
        adjointstateexplicit=retrieveimpulsemodelinformation(ocStruct,'adjointstateexplicit',arcidentifier,symkernel);
        if isempty(adjointstateexplicit.value)
            data.value=[];
            return
        end
        statedynamicsexplicit=retrieveimpulsemodelinformation(ocStruct,'statedynamics');
        optimalvalue=retrieveimpulsemodelinformation(ocStruct,'optimalvalue',arcidentifier);
        maximizingvariable=retrieveimpulsemodelinformation(ocStruct,'maximizingvariable');
        costatenametp1=retrieveimpulsemodelinformation(ocStruct,'costatenametp1');
        statenum=retrieveimpulsemodelinformation(ocStruct,'statenum');
        for ii=1:numel(maximizingvariable.value)
            for jj=1:numel(statedynamicsexplicit.value)
                for kk=1:statenum.value
                    statedynamicsexplicit.value{jj}=ocmatsimple(ocmatsubs(ocmatsubs(statedynamicsexplicit.value{jj},[maximizingvariable.value{ii} '=' optimalvalue.value.(maximizingvariable.value{ii})],symkernel),[costatenametp1.value{kk} '=' adjointstateexplicit.value{kk}],symkernel),symkernel);
                end
            end
        end
        data.value=[statedynamicsexplicit.value adjointstateexplicit.value];
        data.type='mathcellchar';
        data.description='state dynamics, adjoint system (and algebraic equations for implicit control values)';

    case 'specificcanonicalsystemimplicit'
        % the canonical system consisting of the state dynamics, adjoint
        % system and possible algebraic equations for implicitly given
        % control values.
        canonicalsystemimplicit=retrieveimpulsemodelinformation(ocStruct,'canonicalsystemimplicit','',symkernel);
        optimalvalue=retrieveimpulsemodelinformation(ocStruct,'optimalvalue',arcidentifier);
        maximizingvariable=retrieveimpulsemodelinformation(ocStruct,'maximizingvariable');
        canonicalsystem=canonicalsystemimplicit.value;
        for jj=1:numel(canonicalsystem)
            for ii=1:numel(maximizingvariable.value)
                canonicalsystem{jj}=ocmatsimple(ocmatsubs(canonicalsystem{jj},[maximizingvariable.value{ii} '=' optimalvalue.value.(maximizingvariable.value{ii})],symkernel),symkernel);
            end
        end
        data.value=canonicalsystem;
        data.type='mathcellchar';
        data.description='state dynamics, adjoint system (and algebraic equations for implicit control values)';

    case 'specificpontryaginfunction'
        % the canonical system consisting of the state dynamics, adjoint
        % system and possible algebraic equations for implicitly given
        % control values.
        pontryaginfunction=retrieveimpulsemodelinformation(ocStruct,'pontryaginfunction','',symkernel);
        optimalvalue=retrieveimpulsemodelinformation(ocStruct,'optimalvalue',arcidentifier,symkernel);
        maximizingvariable=retrieveimpulsemodelinformation(ocStruct,'maximizingvariable');
        P=pontryaginfunction.value;
        for ii=1:numel(maximizingvariable.value)
            P=ocmatsubs(P,[maximizingvariable.value{ii} '=' optimalvalue.value.(maximizingvariable.value{ii})]);
        end
        data.value=P;
        data.type='mathcellchar';
        data.description='state dynamics, adjoint system (and algebraic equations for implicit control values)';

    case 'canonicalsystemexplicit'
        % the canonical system consisting of the state dynamics, adjoint
        % system and possible algebraic equations for implicitly given
        % control values.
        adjointsystem=retrieveimpulsemodelinformation(ocStruct,'adjointsystemexplicit',arcidentifier,symkernel);
        statedynamics=retrieveimpulsemodelinformation(ocStruct,'statedynamics');
        xp1=[statedynamics.value adjointsystem.value];
        data.value=xp1;
        data.type='mathcellchar';
        data.description='state dynamics, adjoint system';

    case 'canonicalsystemimplicit'
        % the canonical system consisting of the state dynamics, adjoint
        % system and possible algebraic equations for implicitly given
        % control values.
        adjointsystem=retrieveimpulsemodelinformation(ocStruct,'adjointsystem','',symkernel);
        statedynamicsimplicit=retrieveimpulsemodelinformation(ocStruct,'statedynamicsimplicit');
        xp1=[statedynamicsimplicit.value adjointsystem.value];
        data.value=xp1;
        data.type='mathcellchar';
        data.description='state dynamics, adjoint system';

    case 'canonicalsystem'
        % the canonical system consisting of the state dynamics, adjoint
        % system and possible algebraic equations for implicitly given
        % control values.
        adjointsystem=retrieveimpulsemodelinformation(ocStruct,'adjointsystem','',symkernel);
        statedynamics=retrieveimpulsemodelinformation(ocStruct,'statedynamics');
        dxdt=[statedynamics.value adjointsystem.value];
        data.value=dxdt;
        data.type='mathcellchar';
        data.description='state dynamics, adjoint system';

    case 'generalcanonicalsystem'
        % the canonical system consisting of the state dynamics, adjoint
        % system and possible algebraic equations for implicitly given
        % control values.
        adjointsystem=retrieveimpulsemodelinformation(ocStruct,'adjointsystem','',symkernel);
        algebraicequation=retrieveimpulsemodelinformation(ocStruct,'algebraicequation',arcidentifier,symkernel);
        statedynamics=retrieveimpulsemodelinformation(ocStruct,'statedynamics');
        dxdt=[statedynamics.value adjointsystem.value algebraicequation.value];
        data.value=dxdt;
        data.type='mathcellchar';
        data.description='state dynamics, adjoint system (and algebraic equations for implicit control values)';

    case 'canonicalsystemjacobian'
        % Jacobian of the canonical system with respect to state, costate
        % and implicit control/Lagrange multiplier variables
        data.value=ocStruct.foc.canonicalsystem.derivative.(arcidentifier2field(arcidentifier)).DX.term;
        data.type='mathcellchar';
        data.description='Jacobian of the canonical system with respect to state, costate and implicit control/Lagrange multiplier variables';

    case 'canonicalsystemjacobianstatecontrol'
        % Jacobian of the canonical system with respect to state, costate
        % and implicit control/Lagrange multiplier variables
        canonicalsystem=retrieveimpulsemodelinformation(ocStruct,'specificcanonicalsysteminstatecontrol',arcidentifier,symkernel);
        dxdt=cell2vectorstring(canonicalsystem.value);
        equationvariablename=retrieveimpulsemodelinformation(ocStruct,'equationvariablename',arcidentifier);
        variable=cell2vectorstring(equationvariablename.value);
        data.value=string2cell(removematrixstring(ocmatjacobian(dxdt,variable,symkernel)),'matrix');
        data.type='mathcellchar';
        data.description='Jacobian of the canonical system with respect to state, costate and control/Lagrange multiplier variables, with explicit control dynamics';

    case 'statecostatejacobian'
        % Jacobian of the canonical system with respect to state, costate
        % and implicit control/Lagrange multiplier variables
        specificstatedynamics=retrieveimpulsemodelinformation(ocStruct,'specificstatedynamics',arcidentifier,symkernel);
        specificadjointsystem=retrieveimpulsemodelinformation(ocStruct,'specificadjointsystem',arcidentifier,symkernel);
        dxdt=cell2vectorstring([specificstatedynamics.value specificadjointsystem.value]);
        equationvariablename=retrieveimpulsemodelinformation(ocStruct,'equationvariablename',arcidentifier);
        variable=cell2vectorstring(equationvariablename.value);
        data.value=string2cell(removematrixstring(ocmatjacobian(dxdt,variable,symkernel)),'matrix');
        data.type='mathcellchar';
        data.description='Jacobian of the canonical system with respect to state, costate and implicit control/Lagrange multiplier variables';

    case 'canonicalsystemhessian'
        % Jacobian of the canonical system with respect to state, costate
        % and implicit control/Lagrange multiplier variables
        canonicalsystem=retrieveimpulsemodelinformation(ocStruct,'specificcanonicalsystem',arcidentifier,symkernel);
        equationvariablename=retrieveimpulsemodelinformation(ocStruct,'equationvariablename',arcidentifier);
        variable=cell2vectorstring(equationvariablename.value);
        dxdt=canonicalsystem.value;
        for ii=1:numel(dxdt)
            H{ii}=string2cell(removematrixstring(ocmathessian(dxdt{ii},variable,symkernel)),'matrix');
        end
        sx=numel(H);
        if sx
            sy=numel(H{1});
        else
            sy=0;
        end
        for ii=1:sx
            for jj=1:sy
                data.value{jj}{ii}=H{ii}{jj};
            end
        end
        data.type='mathcellchar';
        data.description='Hessian of the canonical system with respect to state, costate and implicit control/Lagrange multiplier variables';

    case 'canonicalsystemparameterjacobian'
        % Jacobian of the canonical system with respect to exogenous
        % parameter variables
        data.value=ocStruct.foc.canonicalsystem.derivative.(arcidentifier2field(arcidentifier)).DP.term;
        data.type='mathcellchar';
        data.description='Jacobian of the canonical system with respect to exogenous parameter variables';
    case 'statecontrolcanonicalsystemparameterjacobian'
        % Jacobian of the canonical system with respect to exogenous
        % parameter variables
        canonicalsystem=retrieveimpulsemodelinformation(ocStruct,'specificcanonicalsysteminstatecontrol',arcidentifier,symkernel);
        dxdt=cell2vectorstring(canonicalsystem.value);
        parametername=retrieveimpulsemodelinformation(ocStruct,'parametername',arcidentifier);
        variable=cell2vectorstring(parametername.value);
        data.value=string2cell(removematrixstring(ocmatjacobian(dxdt,variable,symkernel)),'matrix');
        data.type='mathcellchar';
        data.description='Jacobian of the canonical system with respect to exogenous parameter variables';

    case 'canonicalsystemderivativetime'
        % Jacobian of the canonical system with respect to exogenous
        % parameter variables
        canonicalsystem=retrieveimpulsemodelinformation(ocStruct,'specificcanonicalsystem',arcidentifier,symkernel);
        dxdt=cell2vectorstring(canonicalsystem.value);
        timevariable=retrieveimpulsemodelinformation(ocStruct,'independent');
        variable=cell2vectorstring(timevariable.value);
        data.value=string2cell(removematrixstring(ocmatjacobian(dxdt,variable,symkernel)),'matrix');
        data.type='mathcellchar';
        data.description='Derivative of the canonical system with respect to the independent variable';

    case 'canonicalsystemtotalhessian'
        % Jacobian of the canonical system with respect to exogenous
        % parameter variables
        canonicalsystem=retrieveimpulsemodelinformation(ocStruct,'specificcanonicalsystem',arcidentifier,symkernel);
        equationvariablename=retrieveimpulsemodelinformation(ocStruct,'equationvariablename',arcidentifier);
        parametername=retrieveimpulsemodelinformation(ocStruct,'parametername',arcidentifier);
        variable=cell2vectorstring([equationvariablename.value(:);parametername.value]);
        dxdt=canonicalsystem.value;
        for ii=1:numel(dxdt)
            H{ii}=string2cell(removematrixstring(ocmathessian(dxdt{ii},variable,symkernel)),'matrix');
        end
        sx=numel(H);
        if sx
            sy=numel(H{1});
        else
            sy=0;
        end
        for ii=1:sx
            for jj=1:sy
                data.value{jj}{ii}=H{ii}{jj};
            end
        end
        data.type='mathcellchar';
        data.description='Hessian of the canonical system with respect to state, costate and implicit control/Lagrange multiplier variables';

    case 'optimalvalue'
        % the solutions from the Hamiltonian maximizing condition have
        % already to be stored in the ocStruct structure.
        controlname=retrieveimpulsemodelinformation(ocStruct,'controlname');
        controlnum=retrieveimpulsemodelinformation(ocStruct,'controlnum');
        if controlnum.value
            lagrangemultipliercontrolname=retrieveimpulsemodelinformation(ocStruct,'lagrangemultipliercontrolname');
            inequalitycontrolconstraintnum=retrieveimpulsemodelinformation(ocStruct,'inequalitycontrolconstraintnum');

            for ii=1:controlnum.value
                value.(controlname.value{ii})=ocStruct.foc.value.control.(arcidentifier2field(arcidentifier)).term{ii};
            end
            for ii=1:inequalitycontrolconstraintnum.value
                value.(lagrangemultipliercontrolname.value{ii})=ocStruct.foc.value.lagrangemultcc.(arcidentifier2field(arcidentifier)).term{ii};
            end
        else
            value=[];
        end
        data.type='structmathchar';
        data.value=value;
        data.description='optimal control values as derived from the Hamiltonian maximizing condition';
    case 'impulseoptimalvalue'
        % the solutions from the Hamiltonian maximizing condition have
        % already to be stored in the ocStruct structure.
        icontrolname=retrieveimpulsemodelinformation(ocStruct,'icontrolname');
        icontrolnum=retrieveimpulsemodelinformation(ocStruct,'icontrolnum');
        if icontrolnum.value
            lagrangemultipliercontrolname=retrieveimpulsemodelinformation(ocStruct,'lagrangemultipliercontrolname');
            inequalitycontrolconstraintnum=retrieveimpulsemodelinformation(ocStruct,'inequalitycontrolconstraintnum');

            for ii=1:icontrolnum.value
                value.(icontrolname.value{ii})=ocStruct.foc.value.icontrol.(arcidentifier2field(arcidentifier)).term{ii};
            end
            for ii=1:inequalitycontrolconstraintnum.value
                value.(lagrangemultipliercontrolname.value{ii})=ocStruct.foc.value.lagrangemultcc.(arcidentifier2field(arcidentifier)).term{ii};
            end
        else
            value=[];
        end
        data.type='structmathchar';
        data.value=value;
        data.description='optimal impulse control values as derived from the Hamiltonian maximizing condition';
    case 'hamiltonianDu'
        % derivative of the Pontryaginfunction with respect to the control
        % variable(s)
        controlname=retrieveimpulsemodelinformation(ocStruct,'controlname');
        controlvector=cell2vectorstring(controlname.value);
        hamiltonianfunction=retrieveimpulsemodelinformation(ocStruct,'hamiltonianfunction');
        data.value=string2cell(ocmatjacobian(['[' hamiltonianfunction.value ']'],controlvector,symkernel),'vector');
        data.type='mathcellchar';
        data.description='derivative of the Pontryaginfunction with respect to the control';
    case 'transversalitycondition'
        % derivative of the Pontryaginfunction with respect to the state
        % variable(s)
        %         if isfield(ocStruct,'pontryaginfunction') && isfield(ocStruct.pontryaginfunction,'derivative') ...
        %                 && isfield(ocStruct.pontryaginfunction.derivative,'Dx')
        %             data.value=ocStruct.pontryaginfunction.derivative.Dx.term;
        %         else
        statename=retrieveimpulsemodelinformation(ocStruct,'statenamer');
        statevector=cell2vectorstring(statename.value);
        salvagevalue=retrieveimpulsemodelinformation(ocStruct,'salvagevalue');
        if isempty(salvagevalue.value)
            salvagevalue.value='0';
        end
        data.value=string2cell(ocmatjacobian(['[' salvagevalue.value ']'],statevector,symkernel),'vector');
        %        end
        data.type='mathcellchar';
        data.description='derivative of the salvage value with respect to the state';
    otherwise
        ocmatmsg('No value returned. Property class ''%s'' is unknown.\n',propertyclass)
end


function cellstring=string2cell(string,matrixtype)
% each cell for matrixtype 'matrix' consists of a row of the matrix
% each cell for matrixtype 'vector' consists of an entry of the vector

cellstring=[];
string=removematrixstring(string);
if ~strcmp(string([1 end]),'[]')
    cellstring{1}=string;
    return
end
switch matrixtype
    case 'vector'
        cellstring=regexp(string(2:end-1),',','split');
    case 'matrix'
        cellstring=regexp(string(2:end-1),'],(\ )*[','split');
    case 'charmatrix'
        cellstring=regexp(string(2:end-1),'],[','split');
end

function symval=mycell2sym(cellvec,matrixtype)
% each cell for matrixtype 'matrix' consists of a row of the matrix
% each cell for matrixtype 'vector' consists of an entry of the vector

symval=sym([]);
l=length(cellvec);
switch matrixtype
    case 'vector'
        stringval=['[' deblank(cellvec{1})];
        for ii=2:l
            stringval=[stringval ';' deblank(cellvec{ii})];
        end
        symval=sym([stringval ']']);
    case 'matrix'
        stringval=['[' deblank(cellvec{1})];
        for ii=2:l
            stringval=[stringval ';' deblank(cellvec{ii})];
        end
        symval=sym([stringval ']']);
end

