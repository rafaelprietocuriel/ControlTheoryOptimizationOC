function data=retrievemodelinformation(ocStruct,propertyclass,arcidentifier,symkernel)
%
% RETRIEVEMODELINFORMATION returns a structure with information about the
% ocmat model
% this function is the interface to the structure OCSTRUCT as derived
% from the models initialization file. Its purpose is to allow the same
% commands even if the structure OCSTRUCT is changed. In that case only
% 'retrievemodelinformation' has to be adapted.

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
        constraintcombination=retrievemodelinformation(ocStruct,'constraintcombination',arcidentifier);
        inequalitycontrolconstraintidentifier=retrievemodelinformation(ocStruct,'inequalitycontrolconstraintidentifier');
        inequalitycontrolconstraintnum=retrievemodelinformation(ocStruct,'inequalitycontrolconstraintnum');
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

    case 'modelname'
        % variable name of the state(s) as defined by the user in the
        % initialization file
        data.type='char';
        data.value=ocStruct.modelname;
        data.description='name of the model';

    case 'statename'
        % variable name of the state(s) as defined by the user in the
        % initialization file
        data.type='cellchar';
        data.value=ocStruct.variable.state.name;
        data.description='name of state';

    case 'costatename'
        % variable name of the costate(s) as defined by the user in the
        % initialization file
        data.type='cellchar';
        data.value=ocStruct.variable.costate.name;
        data.description='name of costate';

    case 'controlname'
        % variable name of the control(s) as defined by the user in the
        % initialization file
        data.type='cellchar';
        data.value=ocStruct.variable.control.name;
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
            exogenousfunctionname=retrievemodelinformation(ocStruct,'exogenousfunctionname');
            exogenousfunctionargument=retrievemodelinformation(ocStruct,'exogenousfunctionargument');
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
        exogenousfunctionname=retrievemodelinformation(ocStruct,'exogenousfunctionname');
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
        exogenousfunctionname=retrievemodelinformation(ocStruct,'exogenousfunctionname');
        exogenousfunctionargument=retrievemodelinformation(ocStruct,'exogenousfunctionargument');
        statename=retrievemodelinformation(ocStruct,'statename');
        statenum=retrievemodelinformation(ocStruct,'statenum');
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
        exogenousfunctionname=retrievemodelinformation(ocStruct,'exogenousfunctionname');
        exogenousfunctionargument=retrievemodelinformation(ocStruct,'exogenousfunctionargument');
        statename=retrievemodelinformation(ocStruct,'statename');
        statenum=retrievemodelinformation(ocStruct,'statenum');
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
        exogenousfunctionname=retrievemodelinformation(ocStruct,'exogenousfunctionname');
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
        parametername=retrievemodelinformation(ocStruct,'parametername');
        parametervaluetermindex=retrievemodelinformation(ocStruct,'parametervaluetermindex');
        parameteralgebraictermindex=retrievemodelinformation(ocStruct,'parameteralgebraictermindex');
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

    case 'parametervalueexpression'
        % vector of user provided parameter values
        data.type='cellchar';
        parametername=retrievemodelinformation(ocStruct,'parametername');
        for ii=1:numel(parametername.value)
            if ~isnumeric(ocStruct.parameter.variable.(parametername.value{ii}))
                data.value{ii}=ocStruct.parameter.variable.(parametername.value{ii});
            else
                data.value{ii}=[];
            end
        end
        data.description='vector of user provided parameter value terms';

    case 'discountrate'
        data.type='double';
        parametername=retrievemodelinformation(ocStruct,'parametername');
        parametervalue=retrievemodelinformation(ocStruct,'parametervalue');
        discountratevariable=retrievemodelinformation(ocStruct,'discountratevariable');
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

    case 'lagrangemultiplierstatename'
        % variable name of the Lagrange multiplier(s) for inequality
        % (control) constraints as defined by the user in the
        % initialization file
        if isfield(ocStruct.variable,'lagrangemultsc')
            data.type='cellchar';
            data.value=ocStruct.variable.lagrangemultsc.name;
        end
        data.description='name of Lagrange multiplier for state control';

    case 'inequalitycontrolconstraintidentifier'
        if isfield(ocStruct.constraint,'function') && ...
                isfield(ocStruct.constraint.function,'control') && ...
                ocStruct.constraint.function.control.num
            data.value=ocStruct.constraint.function.control.identifier;
            data.type='cellchar';
        end

    case 'inequalitystateconstraintidentifier'
        if isfield(ocStruct.constraint,'function') && ...
                isfield(ocStruct.constraint.function,'state') && ...
                ocStruct.constraint.function.state.num
            data.value=ocStruct.constraint.function.state.identifier;
            data.type='cellchar';
        end

    case 'arcwithinequalitystateconstraint'
        data.type='char';
        data.value='';
        for arc=ocStruct.arc.identifier
            nonzerolmscindex=retrievemodelinformation(ocStruct,'nonzerolmscindex',arc);
            if ~isempty(nonzerolmscindex.value)
                data.value=[data.value arc{1} ','];
            end
        end
        if ~isempty(data.value)
            data.value(end)=[];
        end
        data.description='number of state';
        
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

    case 'controlnum'
        % the number of control(s)
        data.type='integer';
        data.value=ocStruct.variable.control.num;
        data.description='number of control';

    case 'arcnum'
        % the number of different arcs, corresponding to different
        % functional specifications of the canonical system
        data.type='integer';
        data.value=ocStruct.arc.num;
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
        parameternum=retrievemodelinformation(ocStruct,'parameternum');
        parameteralgebraictermindex=retrievemodelinformation(ocStruct,'parameteralgebraictermindex');
        if isfield(ocStruct,'parameter')
            data.value=setdiff(1:parameternum.value,parameteralgebraictermindex.value);
        else
            data.value=[];
        end
        data.description='index of parameters given by an explicit value';

    case 'odedim'
        % the number of state(s)
        data.type='integer';
        data.value=ocStruct.variable.state.num+ocStruct.variable.costate.num;
        data.description='number of ODEs for the canonical system';

    case 'odestateindex'
        % the number of state(s)
        data.type='integer';
        try
            data.value=find(strcmp(ocStruct.constraint.derivative.state.type,'ode'));
        catch
            data.value=1:ocStruct.variable.state.num;
        end
        data.description='number of ODEs for the canonical system';

    case 'integralstateindex'
        % the number of state(s)
        data.type='integer';
        try
            data.value=find(strcmp(ocStruct.constraint.derivative.state.type,'int'));
        catch
            data.value=0;
        end
        data.description='number of ODEs for the canonical system';
        
    case 'localodestateindex'
        % the number of state(s)
        data.type='integer';
        localstateindex=retrievemodelinformation(ocStruct,'localstateindex');
        odestateindex=retrievemodelinformation(ocStruct,'odestateindex');
        data.value=intersect(localstateindex.value,odestateindex.value);
        data.description='number of ODEs for the canonical system';

    case 'localodestatenum'
        % the number of state(s)
        data.type='integer';
        localodestateindex=retrievemodelinformation(ocStruct,'localodestateindex');
        data.value=numel(localodestateindex.value);
        data.description='number of local states';
        
    case 'nonlocalodestateindex'
        % the number of state(s)
        data.type='integer';
        nonlocalstateindex=retrievemodelinformation(ocStruct,'nonlocalstateindex');
        odestateindex=retrievemodelinformation(ocStruct,'odestateindex');
        data.value=intersect(nonlocalstateindex.value,odestateindex.value);
        data.description='number of ODEs for the canonical system';

    case 'nonlocalodestatenum'
        % the number of state(s)
        data.type='integer';
        nonlocalodestateindex=retrievemodelinformation(ocStruct,'nonlocalodestateindex');
        data.value=numel(nonlocalodestateindex.value);
        data.description='number of local states';

        
    case 'localintegralstateindex'
        % the number of state(s)
        data.type='integer';
        localstateindex=retrievemodelinformation(ocStruct,'localstateindex');
        integralstateindex=retrievemodelinformation(ocStruct,'integralstateindex');
        data.value=intersect(localstateindex.value,integralstateindex.value);
        data.description='number of ODEs for the canonical system';

    case 'localintegralstatenum'
        % the number of state(s)
        data.type='integer';
        localintegralstateindex=retrievemodelinformation(ocStruct,'localintegralstateindex');
        data.value=numel(localintegralstateindex.value);
        data.description='number of local states';
        
    case 'nonlocalintegralstateindex'
        % the number of state(s)
        data.type='integer';
        nonlocalstateindex=retrievemodelinformation(ocStruct,'nonlocalstateindex');
        integralstateindex=retrievemodelinformation(ocStruct,'integralstateindex');
        data.value=intersect(nonlocalstateindex.value,integralstateindex.value);
        data.description='number of ODEs for the canonical system';

    case 'nonlocalintegralstatenum'
        % the number of state(s)
        data.type='integer';
        nonlocalintegralstateindex=retrievemodelinformation(ocStruct,'nonlocalintegralstateindex');
        data.value=numel(nonlocalintegralstateindex.value);
        data.description='number of local states';
        
    case 'localstatenum'
        % the number of state(s)
        data.type='integer';
        localstatenum=retrievemodelinformation(ocStruct,'localstatenum');
        data.value=numel(localstatenum.value);
        data.description='number of local states';

    case 'localstateindex'
        % the number of state(s)
        data.type='integer';
        try
            data.value=find(strcmp(ocStruct.variable.state.property,'l'));
        catch
            data.value=1:ocStruct.variable.state.num;
        end
        data.description='index of local states';

    case 'nonlocalstatenum'
        % the number of state(s)
        data.type='integer';
        nonlocalstatenum=retrievemodelinformation(ocStruct,'nonlocalstatenum');
        data.value=numel(nonlocalstatenum.value);
        data.description='number of state';

    case 'nonlocalstateindex'
        % the number of state(s)
        data.type='integer';
        try
            data.value=find(strcmp(ocStruct.variable.state.property,'nl'));
        catch
            data.value=1:ocStruct.variable.state.num;
        end
        data.description='index of non-local states';

    case 'localcontrolnum'
        % the number of state(s)
        data.type='integer';
        localcontrolindex=retrievemodelinformation(ocStruct,'localcontrolindex');
        data.value=numel(localcontrolindex.value);
        data.description='number of local controls';

    case 'localcontrolindex'
        % the number of control(s)
        data.type='integer';
        try
            data.value=find(strcmp(ocStruct.variable.control.property,'l'));
        catch
            data.value=0;
        end
        data.description='index of local controls';

    case 'nonlocalcontrolnum'
        % the number of control(s)
        data.type='integer';
        nonlocalcontrolindex=retrievemodelinformation(ocStruct,'nonlocalcontrolindex');
        data.value=numel(nonlocalcontrolindex.value);
        data.description='number of control';

    case 'nonlocalcontrolindex'
        % the number of control(s)
        data.type='integer';
        try
            data.value=find(strcmp(ocStruct.variable.control.property,'nl'));
        catch
            data.value=1:ocStruct.variable.control.num;
        end
        data.description='index of non-local controls';
        
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
            data.value=length(ocStruct.constraint.function.control);
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

    case 'inequalitystateconstrainttimederivative'
        % inequality constraints which only include state variables
        if isfield(ocStruct.constraint,'function') &&  ...
                isfield(ocStruct.constraint.function,'state') && ...
                ocStruct.constraint.function.state.num

            data.value=ocStruct.constraint.function.state.timederivative;
            data.type='mathchar';
        end
        data.description='time derivative pure state inequality constraints';

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

    case 'inequalitystateconstraintorder'
        inequalitystateconstraintnum=retrievemodelinformation(ocStruct,'inequalitystateconstraintnum');
        data.type='integer vector';
        if inequalitystateconstraintnum.value
            data.value=[ocStruct.constraint.function.state.order{:}];
        else
            data.value=[];
        end
        data.description='order of pure state inequality constraints';

    case 'algebraicequationnum'
        maximizingimplicitderivativevariable=retrievemodelinformation(ocStruct,'maximizingimplicitderivativevariable',arcidentifier);
        data.type='integer';
        data.value=numel(maximizingimplicitderivativevariable.value);
        data.description='number of algebraic equations';

    case 'totalalgebraicequationnum'
        totalimplicitvariable=retrievemodelinformation(ocStruct,'totalimplicitvariable');
        data.type='integer';
        data.value=numel(totalimplicitvariable.value);
        data.description='number of algebraic equations';

    case 'canonicalsystemequationnum'
        statenum=retrievemodelinformation(ocStruct,'statenum');
        costatenum=retrievemodelinformation(ocStruct,'costatenum');
        algebraicequationnum=retrievemodelinformation(ocStruct,'algebraicequationnum',arcidentifier);
        data.type='integer';
        data.value=statenum.value+costatenum.value+algebraicequationnum.value;
        data.description='total number of algebraic differential equations';

    case 'statedynamics'
        % ODEs describing the evolution of the state(s)
        data.type='mathchar';
        data.value=ocStruct.constraint.derivative.state.term;
        data.description='ODEs of the state dynamics';

    case 'objectivetype'
        % type of the terms which build up the objective value function
        % (usually integral and Salvage value)
        data.type='mathchar';
        data.value=fieldnames(ocStruct.objective);
        data.description='type of terms which build up the objective value function';

    case 'objectiveintegrand'
        % integrand for objective value function without the possible
        % discounting factor
        data.type='mathchar';
        data.value=ocStruct.objective.integral.function.term;
        data.description='integrand for objective value function';

    case 'salvagevalue'
        % Salvagevalue
        data.type='mathchar';
        if isfield(ocStruct.objective,'sum') && isfield(ocStruct.objective.sum,'endtime')
            data.value=ocStruct.objective.sum.endtime.function.term;
        else
            data.value='0';
        end
        data.description='salvagevalue';

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

    case 'discobjectivefunction'
        objectivefunction=retrievemodelinformation(ocStruct,'objectiveintegrand');
        discountratevariable=retrievemodelinformation(ocStruct,'discountratevariable');
        independentvariable=retrievemodelinformation(ocStruct,'independent');
        if isfield(ocStruct.objective,'integral')
            data.type='mathchar';
            data.value=['exp(-' discountratevariable.value '*' independentvariable.value ')*(' objectivefunction.value ')'];
            data.description='the discounted objective integrand';
        end

    case 'discobjectivefunctionderivativetime'
        specificdiscobjectivefunction=retrievemodelinformation(ocStruct,'specificdiscobjectivefunction',arcidentifier,symkernel);
        timevariable=retrievemodelinformation(ocStruct,'independent');
        variable=cell2vectorstring(timevariable.value);
        data.value=string2cell(removematrixstring(ocmatdiff(specificdiscobjectivefunction.value,variable,symkernel)));
        data.type='mathchar';
        data.description='the discounted objective integrand';

    case 'specificdiscobjectivefunction'
        discobjectivefunction=retrievemodelinformation(ocStruct,'discobjectivefunction');
        optimalvalue=retrievemodelinformation(ocStruct,'optimalvalue',arcidentifier);
        maximizingvariable=retrievemodelinformation(ocStruct,'maximizingvariable');
        for ii=1:numel(maximizingvariable.value)
            discobjectivefunction.value=ocmatsubs(discobjectivefunction.value,[maximizingvariable.value{ii} '=' optimalvalue.value.(maximizingvariable.value{ii})],symkernel);
        end
        data.type='mathchar';
        data.value=discobjectivefunction.value;
        data.description='the discounted objective integrand';

    case 'discobjectivefunctionDX'
        costatename=retrievemodelinformation(ocStruct,'costatename');
        statename=retrievemodelinformation(ocStruct,'statename');
        statevector=cell2vectorstring([statename.value costatename.value]);
        specificdiscobjectivefunction=retrievemodelinformation(ocStruct,'specificdiscobjectivefunction',arcidentifier,symkernel);
        data.type='mathcellchar';
        data.value=string2cell(ocmatjacobian(['[' specificdiscobjectivefunction.value ']'],statevector,symkernel),'vector');
        data.description='derivative of the discounted objective integrand with respect to the state and costate';

    case 'discobjectivefunctionDP'
        specificdiscobjectivefunction=retrievemodelinformation(ocStruct,'specificdiscobjectivefunction',arcidentifier,symkernel);
        parametername=retrievemodelinformation(ocStruct,'parametername',arcidentifier);
        variable=cell2vectorstring(parametername.value);
        data.type='mathcellchar';
        data.value=string2cell(ocmatjacobian(['[' specificdiscobjectivefunction.value ']'],variable,symkernel),'vector');
        data.description='derivative of the discounted objective integrand with respect to the state and costate';

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

    case 'equationvariablename'
        % dependent variables of the canonical system for a specific arc
        data.type='cellchar';
        statename=retrievemodelinformation(ocStruct,'statename');
        costatename=retrievemodelinformation(ocStruct,'costatename');
        implicitcontrolname=retrievemodelinformation(ocStruct,'implicitnonlinearcontrol',arcidentifier);
        equationvariablename=[statename.value costatename.value implicitcontrolname.value];
        transform2statecontrolspace=retrievemodelinformation(ocStruct,'transform2statecontrolspace',arcidentifier);
        if transform2statecontrolspace.value
            replacedcostatename=retrievemodelinformation(ocStruct,'replacedcostatename',arcidentifier);
            controlname4controlstate=retrievemodelinformation(ocStruct,'controlname4controlstate',arcidentifier);
            for ii=1:numel(replacedcostatename.value)
                equationvariablename=strrep(equationvariablename,replacedcostatename.value{ii},controlname4controlstate.value{ii});
            end
        end
        data.value=equationvariablename;
        data.description='dependent variables of the canonical system';

    case 'explicitnonlinearcontrol'
        % explicit nonlinear control variables for a specific arc
        property=ocStruct.variable.control.(arcidentifier2field(arcidentifier)).property;
        data.value=ocStruct.variable.control.name([property.explicit{:}]==1 & [property.linear{:}]==0);
        data.type='cellchar';
        data.description='explicit nonlinear control variable name';

    case 'explicitnonlinearcontrolindex'
        % explicit nonlinear control variables for a specific arc
        property=ocStruct.variable.control.(arcidentifier2field(arcidentifier)).property;
        data.value=find([property.explicit{:}]==1 & [property.linear{:}]==0);
        data.type='cellchar';
        data.description='explicit nonlinear control variable name';

    case 'implicitnonlinearcontrol'
        % implicit nonlinear control variables for a specific arc
        property=ocStruct.variable.control.(arcidentifier2field(arcidentifier)).property;
        data.value=ocStruct.variable.control.name([property.explicit{:}]==0 & [property.linear{:}]==0);
        data.type='cellchar';
        data.description='implicit nonlinear control variable name';

    case 'implicitnonlinearcontrolindex'
        % implicit nonlinear control variables for a specific arc
        property=ocStruct.variable.control.(arcidentifier2field(arcidentifier)).property;
        data.value=find([property.explicit{:}]==0 & [property.linear{:}]==0);
        data.type='integer';
        data.description='implicit nonlinear control variable indices';

    case 'transform2statecontrolspace'
        % implicit nonlinear control variables for a specific arc
        property=ocStruct.variable.control.(arcidentifier2field(arcidentifier)).property;
        if isfield(property,'replace')
            data.value=any([property.replace{:}]~=0);
        else
            data.value=0;
        end
        data.type='integer';
        data.description='index of control variable for which the control dynamics has to be derived';

    case 'explicitcontroldynamicsindex'
        % implicit nonlinear control variables for a specific arc
        property=ocStruct.variable.control.(arcidentifier2field(arcidentifier)).property;
        if isfield(property,'replace')
            data.value=find([property.replace{:}]~=0);
        else
            data.value=[];
        end
        data.type='integer';
        data.description='index of control variable for which the control dynamics has to be derived';

    case 'replacedcostateindex'
        % implicit nonlinear control variables for a specific arc
        property=ocStruct.variable.control.(arcidentifier2field(arcidentifier)).property;
        index=[property.replace{:}];
        index(index==0)=[];
        data.value=index;
        data.type='integer';
        data.description='index of control variable for which the control dynamics has to be derived';

    case 'replacedcostatename'
        replacedcostateindex=retrievemodelinformation(ocStruct,'replacedcostateindex',arcidentifier);
        costatename=retrievemodelinformation(ocStruct,'costatename');
        costatename.value=costatename.value(replacedcostateindex.value);
        data.value=costatename.value;
        data.type='cellchar';
        data.description='linear control variables';

    case 'controlname4controlstate'
        explicitcontroldynamicsindex=retrievemodelinformation(ocStruct,'explicitcontroldynamicsindex',arcidentifier);
        controlname=retrievemodelinformation(ocStruct,'controlname');
        controlname.value=controlname.value(explicitcontroldynamicsindex.value);
        data.value=controlname.value;
        data.type='cellchar';
        data.description='linear control variables';

    case 'linearcontrol'
        % linear control variables for a specific arc
        property=ocStruct.variable.control.(arcidentifier2field(arcidentifier)).property;
        data.value=ocStruct.variable.control.name([property.linear{:}]==1);
        data.type='cellchar';
        data.description='linear control variables';

    case 'maximizingvariable'
        % dependent variable names of the Hamiltonian maximizing condition
        controlname=retrievemodelinformation(ocStruct,'controlname');
        lagrangemultipliercontrolname=retrievemodelinformation(ocStruct,'lagrangemultipliercontrolname');
        lagrangemultiplierstatename=retrievemodelinformation(ocStruct,'lagrangemultiplierstatename');
        data.value=[controlname.value lagrangemultipliercontrolname.value lagrangemultiplierstatename.value];
        data.type='cellchar';
        data.description='variable names of the Hamiltonian maximizing condition';

    case 'maximizingderivativevariable'
        % variable names which are used to derive the first order necessary
        % conditions of the Hamiltonian maximizing condition
        arc=arcidentifier2field(arcidentifier);
        if isfield(ocStruct,'foc') && isfield(ocStruct.foc,'generalinformation') && ...
                isfield(ocStruct.foc.generalinformation,arc) && ...
                isfield(ocStruct.foc.generalinformation.(arc),'maximizingderivativevariable')

            data.value=ocStruct.foc.generalinformation.(arc).maximizingderivativevariable.name;
        else
            explicitnonlinearcontrol=retrievemodelinformation(ocStruct,'explicitnonlinearcontrol',arcidentifier);
            nonzerolmmc=retrievemodelinformation(ocStruct,'nonzerolmmc',arcidentifier);
            nonzerolmsc=retrievemodelinformation(ocStruct,'nonzerolmsc',arcidentifier);
            data.value=[explicitnonlinearcontrol.value nonzerolmmc.value nonzerolmsc.value];
        end
        data.type='cellchar';
        data.description='variable names for the first order necessary conditions of the Hamiltonian maximizing condition';

    case 'maximizingimplicitderivativevariable'
        % variable names which are used to derive the first order necessary
        % conditions of the Hamiltonian maximizing condition and are only
        % calculated implicitly
        data=retrievemodelinformation(ocStruct,'maximizingvariable');
        maximizingderivativevariable=retrievemodelinformation(ocStruct,'maximizingderivativevariable',arcidentifier);
        zerolmmc=retrievemodelinformation(ocStruct,'zerolmmc',arcidentifier);
        zerolmsc=retrievemodelinformation(ocStruct,'zerolmsc',arcidentifier);
        totalexplicitvariable=[maximizingderivativevariable.value zerolmmc.value zerolmsc.value];
        for ii=1:numel(totalexplicitvariable)
            data.value(strcmp(data.value,totalexplicitvariable{ii}))=[];
        end
        data.description='implicit variable names for the first order necessary conditions of the Hamiltonian maximizing condition';

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
            explicitnonlinearcontrol=retrievemodelinformation(ocStruct,'explicitnonlinearcontrol',arcidentifier);
            nonzerolmmc=retrievemodelinformation(ocStruct,'nonzerolmmc',arcidentifier);
            nonzerolmsc=retrievemodelinformation(ocStruct,'nonzerolmsc',arcidentifier);
            data.value=[explicitnonlinearcontrol.value nonzerolmmc.value nonzerolmsc.value];
        end
        data.type='cellchar';
        data.description='explicit variable names for the first order necessary conditions of the Hamiltonian maximizing condition';

    case 'maximizingimplicitvariable'
        % variable names which are implicitly calculated by the Hamiltonian
        % maximizing condition
        data=retrievemodelinformation(ocStruct,'maximizingvariable');
        maximizingexplicitvariable=retrievemodelinformation(ocStruct,'maximizingexplicitvariable',arcidentifier);
        zerolmmc=retrievemodelinformation(ocStruct,'zerolmmc',arcidentifier);
        zerolmsc=retrievemodelinformation(ocStruct,'zerolmsc',arcidentifier);
        totalexplicitvariable=[maximizingexplicitvariable.value zerolmmc.value zerolmsc.value];
        for ii=1:numel(totalexplicitvariable)
            data.value(strcmp(data.value,totalexplicitvariable{ii}))=[];
        end
        data.description='implicit variable names of the Hamiltonian maximizing condition';

    case 'totalimplicitvariable'
        controlname=retrievemodelinformation(ocStruct,'controlname');
        totalimplicitnonlinearcontrolindex=retrievemodelinformation(ocStruct,'totalimplicitvariableindex');
        data.value=controlname.value(totalimplicitnonlinearcontrolindex.value);
        data.type='integer';
        data.description='index of the Lagrange multipliers corresponding to an inactive constraint';

    case 'totalimplicitvariableindex'
        arcidentifier=retrievemodelinformation(ocStruct,'arcidentifier');
        totalimplicitnonlinearcontrolindex=[];
        for ii=1:numel(arcidentifier.value)
            implicitnonlinearcontrolindex=retrievemodelinformation(ocStruct,'implicitnonlinearcontrolindex',arcidentifier.value{ii});
            totalimplicitnonlinearcontrolindex=[totalimplicitnonlinearcontrolindex implicitnonlinearcontrolindex.value];
        end
        data.value=unique(totalimplicitnonlinearcontrolindex);
        data.type='integer';
        data.description='index of the Lagrange multipliers corresponding to an inactive constraint';

    case 'equationvariablenameimplicit'
        % dependent variables of the canonical system for a specific arc
        data.type='cellchar';
        statename=retrievemodelinformation(ocStruct,'statename');
        costatename=retrievemodelinformation(ocStruct,'costatename');
        implicitcontrolname=retrievemodelinformation(ocStruct,'maximizingimplicitvariable',arcidentifier);
        data.value=[statename.value costatename.value implicitcontrolname.value];
        data.description='dependent variables of the canonical system';

    case 'zerolmmcindex'
        % returns the index of the Lagrange multipliers corresponding to an
        % inactive constraint for a specific arc
        inequalitycontrolconstraintidentifier=retrievemodelinformation(ocStruct,'inequalitycontrolconstraintidentifier');
        controlconstraintnum=retrievemodelinformation(ocStruct,'inequalitycontrolconstraintnum');
        constraintcombination=retrievemodelinformation(ocStruct,'constraintcombination',arcidentifier);
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
        lagrangemultipliercontrolname=retrievemodelinformation(ocStruct,'lagrangemultipliercontrolname');
        zerolmmcindex=retrievemodelinformation(ocStruct,'zerolmmcindex',arcidentifier);
        data.value=lagrangemultipliercontrolname.value(zerolmmcindex.value);
        data.type='cellchar';
        data.description='name of the Lagrange multipliers corresponding to an inactive constraint';

    case 'nonzerolmmc'
        % returns the name of the Lagrange multipliers corresponding to an
        % active constraint for a specific arc
        lagrangemultipliercontrolname=retrievemodelinformation(ocStruct,'lagrangemultipliercontrolname');
        nonzerolmmcindex=retrievemodelinformation(ocStruct,'nonzerolmmcindex',arcidentifier);
        data.value=lagrangemultipliercontrolname.value(nonzerolmmcindex.value);
        data.type='cellchar';
        data.description='name of the Lagrange multipliers corresponding to an active constraint';

    case 'nonzerolmmcindex'
        % returns the index of the Lagrange multipliers corresponding to an
        % active constraint for a specific arc
        inequalitycontrolconstraintnum=retrievemodelinformation(ocStruct,'inequalitycontrolconstraintnum');
        zerolmmcindex=retrievemodelinformation(ocStruct,'zerolmmcindex',arcidentifier);
        data.value=setdiff(1:inequalitycontrolconstraintnum.value,zerolmmcindex.value);
        data.type='integer';
        data.description='index of the Lagrange multipliers corresponding to an active constraint';


    case 'zerolmscindex'
        % returns the index of the Lagrange multipliers corresponding to an
        % inactive constraint for a specific arc
        inequalitystateconstraintidentifier=retrievemodelinformation(ocStruct,'inequalitystateconstraintidentifier');
        stateconstraintnum=retrievemodelinformation(ocStruct,'inequalitystateconstraintnum');
        constraintcombination=retrievemodelinformation(ocStruct,'constraintcombination',arcidentifier);
        zerolmscindex=[];
        if ~strcmp(constraintcombination.value,'[]')
            for jj=1:numel(inequalitystateconstraintidentifier.value)
                testexpr=regexp(constraintcombination.value,['\<' inequalitystateconstraintidentifier.value{jj} '\>'],'ONCE');
                if isempty([testexpr{:}])
                    zerolmscindex=[zerolmscindex jj];
                end
            end
        else
            zerolmscindex=1:stateconstraintnum.value;
        end
        data.value=zerolmscindex;
        data.type='integer';
        data.description='index of the Lagrange multipliers corresponding to an inactive constraint';

    case 'zerolmsc'
        % returns the index of the Lagrange multipliers corresponding to an
        % inactive constraint for a specific arc
        lagrangemultiplierstatename=retrievemodelinformation(ocStruct,'lagrangemultiplierstatename');
        zerolmscindex=retrievemodelinformation(ocStruct,'zerolmscindex',arcidentifier);
        data.value=lagrangemultiplierstatename.value(zerolmscindex.value);
        data.type='cellchar';
        data.description='name of the Lagrange multipliers corresponding to an inactive constraint';
        
    case 'explicitstatename'
        nonzerolmscindex=retrievemodelinformation(ocStruct,'nonzerolmscindex',arcidentifier);
        if isempty(nonzerolmscindex.value)
            data.value='';
        else
            data.value=ocStruct.variable.state.(arcidentifier2field(arcidentifier)).explicit.name;
        end
        data.type='mathchar';
        data.description='state values determined by the state constraints';

    case 'explicitstateindex'
        nonzerolmscindex=retrievemodelinformation(ocStruct,'nonzerolmscindex',arcidentifier);
        statename=retrievemodelinformation(ocStruct,'statename');
        if isempty(nonzerolmscindex.value)
            data.value=[];
        else
            explicitstatename=retrievemodelinformation(ocStruct,'explicitstatename',arcidentifier);
            for ii=1:numel(explicitstatename.value)
                idx=[];
                tmp=find(strcmp(statename.value,explicitstatename.value{ii}));
                if ~isempty(tmp)
                    idx=[idx tmp];
                end
            end
            data.value=idx;
        end
        data.type='mathchar';
        data.description='index of states explicitly determined by the state constraints';

    case 'explicitstatevalue'
        nonzerolmscindex=retrievemodelinformation(ocStruct,'nonzerolmscindex',arcidentifier);
        statename=retrievemodelinformation(ocStruct,'statename');
        if isempty(nonzerolmscindex.value)
            data.value=statename.value;
        else
            explicitstatename=retrievemodelinformation(ocStruct,'explicitstatename',arcidentifier);
            inequalitystateconstraint=retrievemodelinformation(ocStruct,'inequalitystateconstraint');
            equation=cell2vectorstring([inequalitystateconstraint.value(nonzerolmscindex.value)]);
            equationvariable=cell2vectorstring([explicitstatename.value]);
            solution=ocmatsolve(equation,equationvariable,symkernel);
            if ~isempty(solution)
                for ii=1:numel(solution)
                    % the implicit control values are added to the solution
                    % structure, where the value is the variable itself
                    for jj=1:numel(statename.value)
                        if ~any(strcmp(explicitstatename.value,statename.value{jj}))
                            solution.(statename.value{jj})=statename.value{jj};
                        end
                    end
                end
            end
            solution=orderfields(solution,statename.value);
            for ii=1:length(statename.value)
                data.value{ii}=solution.(statename.value{ii});
            end
        end
        data.type='mathchar';
        data.description='state values determined by the state constraints';
        
    case 'nonzerolmsc'
        % returns the name of the Lagrange multipliers corresponding to an
        % active constraint for a specific arc
        lagrangemultiplierstatename=retrievemodelinformation(ocStruct,'lagrangemultiplierstatename');
        nonzerolmscindex=retrievemodelinformation(ocStruct,'nonzerolmscindex',arcidentifier);
        data.value=lagrangemultiplierstatename.value(nonzerolmscindex.value);
        data.type='cellchar';
        data.description='name of the Lagrange multipliers corresponding to an active constraint';

    case 'nonzerolmscindex'
        % returns the index of the Lagrange multipliers corresponding to an
        % active constraint for a specific arc
        inequalitystateconstraintnum=retrievemodelinformation(ocStruct,'inequalitystateconstraintnum');
        zerolmscindex=retrievemodelinformation(ocStruct,'zerolmscindex',arcidentifier);
        data.value=setdiff(1:inequalitystateconstraintnum.value,zerolmscindex.value);
        data.type='integer';
        data.description='index of the Lagrange multipliers corresponding to an active constraint';

    case 'hamiltonianfunction'
        % general Hamilton function including the "Lagrange terms" from
        % inequality contraints
        if isfield(ocStruct,'hamiltonianfunction') && isfield(ocStruct.hamiltonianfunction,'term')
            data.value=ocStruct.hamiltonianfunction.term;
        else
            statenum=retrievemodelinformation(ocStruct,'statenum');
            costatename=retrievemodelinformation(ocStruct,'costatename');
            % integrand of the objective value function without discount factor
            g=retrievemodelinformation(ocStruct,'objectiveintegrand');
            % state dynamics
            dxdt=retrievemodelinformation(ocStruct,'statedynamics');
            data.value=g.value;
            for ii=1:statenum.value
                data.value=[data.value '+' costatename.value{ii} '*(' dxdt.value{ii} ')'];
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
            statenum=retrievemodelinformation(ocStruct,'statenum');
            costatename=retrievemodelinformation(ocStruct,'costatename');
            inequalitycontrolconstraint=retrievemodelinformation(ocStruct,'inequalitycontrolconstraint');
            inequalitycontrolconstraintnum=retrievemodelinformation(ocStruct,'inequalitycontrolconstraintnum');
            lmmc=retrievemodelinformation(ocStruct,'lagrangemultipliercontrolname');
            inequalitystateconstraintnum=retrievemodelinformation(ocStruct,'inequalitystateconstraintnum');
            if inequalitystateconstraintnum.value
                inequalitystateconstrainttimederivative=retrievemodelinformation(ocStruct,'inequalitystateconstrainttimederivative');
                inequalitystateconstraintorder=retrievemodelinformation(ocStruct,'inequalitystateconstraintorder');
                lmsc=retrievemodelinformation(ocStruct,'lagrangemultiplierstatename');
            end
            % integrand of the objective value function without discount factor
            g=retrievemodelinformation(ocStruct,'objectiveintegrand');
            % state dynamics
            dxdt=retrievemodelinformation(ocStruct,'statedynamics');
            data.value=g.value;
            for ii=1:statenum.value
                data.value=[data.value '+' costatename.value{ii} '*(' dxdt.value{ii} ')'];
            end
            for ii=1:inequalitycontrolconstraintnum.value
                data.value=[data.value '+' lmmc.value{ii} '*(' inequalitycontrolconstraint.value{ii} ')'];
            end
            for ii=1:inequalitystateconstraintnum.value
                data.value=[data.value '+' lmsc.value{ii} '*(' inequalitystateconstrainttimederivative.value{ii}{inequalitystateconstraintorder.value(ii)} ')'];
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
            statename=retrievemodelinformation(ocStruct,'statename');
            statevector=cell2vectorstring(statename.value);
            pontryaginfunction=retrievemodelinformation(ocStruct,'pontryaginfunction');
            data.value=string2cell(ocmatjacobian(['[' pontryaginfunction.value ']'],statevector,symkernel),'vector');
        end
        data.type='mathcellchar';
        data.description='derivative of the Pontryaginfunction with respect to the state';

    case 'pontryaginfunctionDX'
        % derivative of the Pontryaginfunction with respect to the state
        % AND costate variable(s)
        costatename=retrievemodelinformation(ocStruct,'costatename');
        statename=retrievemodelinformation(ocStruct,'statename');
        statevector=cell2vectorstring([statename.value costatename.value]);
        pontryaginfunction=retrievemodelinformation(ocStruct,'specificpontryaginfunction',arcidentifier,symkernel);
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
            controlname=retrievemodelinformation(ocStruct,'controlname');
            controlvector=cell2vectorstring(controlname.value);
            pontryaginfunction=retrievemodelinformation(ocStruct,'pontryaginfunction');
            data.value=string2cell(ocmatjacobian(['[' pontryaginfunction.value ']'],controlvector,symkernel),'vector');
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
            statename=retrievemodelinformation(ocStruct,'statename');
            statevector=cell2vectorstring(statename.value);
            pontryaginfunction=retrievemodelinformation(ocStruct,'pontryaginfunction');
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
            controlname=retrievemodelinformation(ocStruct,'controlname');
            controlvector=cell2vectorstring(controlname.value);
            pontryaginfunction=retrievemodelinformation(ocStruct,'pontryaginfunction');
            data.value=string2cell(ocmathessian(pontryaginfunction.value,controlvector,symkernel),'matrix');
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
            controlname=retrievemodelinformation(ocStruct,'controlname');
            controlvector=cell2vectorstring(controlname.value);
            statename=retrievemodelinformation(ocStruct,'statename');
            statevector=cell2vectorstring(statename.value);
            pontryaginfunction=retrievemodelinformation(ocStruct,'pontryaginfunction');
            J=removematrixstring(ocmatjacobian(['[' pontryaginfunction.value ']'],controlvector,symkernel));
            data.value=string2cell(ocmatjacobian(J,statevector,symkernel),'matrix');
        end
        data.type='mathcellchar';
        data.description='second order derivative of the Pontryaginfunction with respect to the control';

    case 'optimalcontroldynamicsleftside'
        controlname=retrievemodelinformation(ocStruct,'controlname');
        controlvector=cell2vectorstring(controlname.value);
        zerolmmcindex=retrievemodelinformation(ocStruct,'zerolmmcindex',arcidentifier);
        lagrangemultipliercontrolname=retrievemodelinformation(ocStruct,'lagrangemultipliercontrolname');
        lagrangemultipliercontrolname.value(zerolmmcindex.value)=[];
        lagrangemultipliervector=cell2vectorstring(lagrangemultipliercontrolname.value);
        pontryaginfunction=retrievemodelinformation(ocStruct,'pontryaginfunction');
        pontryaginfunction4arc=retrievemodelinformation(ocStruct,'pontryaginfunction4arc',arcidentifier,symkernel);
        pontryaginfunctionDu2=mystr2sym(ocmathessian(pontryaginfunction.value,controlvector,symkernel));
        if ~isempty(lagrangemultipliervector)
            pontryaginfunctionDlmmcDu=mystr2sym(ocmatjacobian(removematrixstring(ocmatjacobian(['[' pontryaginfunction4arc.value ']'],controlvector,symkernel)),lagrangemultipliervector,symkernel));
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

    case 'invoptimalcontroldynamicsleftside'
        optimalcontroldynamicsleftside=retrievemodelinformation(ocStruct,'optimalcontroldynamicsleftside',arcidentifier,symkernel);
        LSM=mycell2sym(optimalcontroldynamicsleftside.value,'matrix');
        data.value=string2cell(char(inv(LSM)),'charmatrix');
        data.type='mathcellchar';
        data.description='inverse of the left side matrix of the optimal control dynamics';


    case 'optimalcontroldynamics'
        invoptimalcontroldynamicsleftside=retrievemodelinformation(ocStruct,'invoptimalcontroldynamicsleftside',arcidentifier,symkernel);
        pontryaginfunctionDuDX=retrievemodelinformation(ocStruct,'pontryaginfunctionDuDX',arcidentifier,symkernel);
        pontryaginfunctionDlmmcDX=retrievemodelinformation(ocStruct,'pontryaginfunctionDlmmcDX',arcidentifier,symkernel);
        %implicitnonlinearcontrolindex=retrievemodelinformation(ocStruct,'implicitnonlinearcontrolindex',arcidentifier);
        implicitnonlinearcontrolindex=retrievemodelinformation(ocStruct,'totalimplicitvariableindex');
        if ~isempty(implicitnonlinearcontrolindex.value)
            canonicalsystem=retrievemodelinformation(ocStruct,'canonicalsystem','',symkernel);
            invLSM=mycell2sym(invoptimalcontroldynamicsleftside.value,'matrix');
            RSM1=mycell2sym(pontryaginfunctionDuDX.value,'matrix');
            if ~isempty(pontryaginfunctionDlmmcDX.value)
                RSM2=mycell2sym(pontryaginfunctionDlmmcDX.value,'matrix');
            else
                RSM2=sym([]);
            end
            RSM=[RSM1;RSM2'];
            dxdt=mycell2sym(canonicalsystem.value,'vector');
            dudt=-invLSM*RSM*dxdt;
            data.value=string2cell(char(dudt(implicitnonlinearcontrolindex.value)),'charmatrix');
            
        else
            data.value='';
        end
        data.type='mathcellchar';
        data.description='inverse of the left side matrix of the optimal control dynamics';

    case 'fulloptimalcontroldynamics'
        invoptimalcontroldynamicsleftside=retrievemodelinformation(ocStruct,'invoptimalcontroldynamicsleftside',arcidentifier,symkernel);
        controlnum=retrievemodelinformation(ocStruct,'controlnum');
        %zerolmmcindex=retrievemodelinformation(ocStruct,'zerolmmcindex',arcidentifier);
        nonzerolmmcindex=retrievemodelinformation(ocStruct,'nonzerolmmcindex',arcidentifier);
        controlconstraintnum=retrievemodelinformation(ocStruct,'inequalitycontrolconstraintnum');

        pontryaginfunctionDuDX=retrievemodelinformation(ocStruct,'pontryaginfunctionDuDX',arcidentifier,symkernel);
        pontryaginfunctionDlmmcDX=retrievemodelinformation(ocStruct,'pontryaginfunctionDlmmcDX',arcidentifier,symkernel);
        canonicalsystem=retrievemodelinformation(ocStruct,'canonicalsystem','',symkernel);
        invLSM=mycell2sym(invoptimalcontroldynamicsleftside.value,'matrix');
        RSM1=mycell2sym(pontryaginfunctionDuDX.value,'matrix');
        if ~isempty(pontryaginfunctionDlmmcDX.value)
            RSM2=mycell2sym(pontryaginfunctionDlmmcDX.value,'matrix');
        else
            RSM2=sym([]);
        end
        RSM=[RSM1;RSM2'];
        dxdt=mycell2sym(canonicalsystem.value,'vector');
        dudt=-invLSM*RSM*dxdt;
        dUdt=sym(zeros(controlnum.value+controlconstraintnum.value,1));
        dUdt([1:controlnum.value controlnum.value+nonzerolmmcindex.value])=dudt;
        data.value=string2cell(char(dUdt),'charmatrix');

        data.type='mathcellchar';
        data.description='dynamics of the controls and Lagrange multipliers';

    case 'explicitoptimalcontroldynamics'
        invoptimalcontroldynamicsleftside=retrievemodelinformation(ocStruct,'invoptimalcontroldynamicsleftside',arcidentifier,symkernel);
        pontryaginfunctionDuDX=retrievemodelinformation(ocStruct,'pontryaginfunctionDuDX',arcidentifier,symkernel);
        pontryaginfunctionDlmmcDX=retrievemodelinformation(ocStruct,'pontryaginfunctionDlmmcDX',arcidentifier,symkernel);
        explicitcontroldynamicsindex=retrievemodelinformation(ocStruct,'explicitcontroldynamicsindex',arcidentifier);
        optimalcostatevalue=retrievemodelinformation(ocStruct,'optimalcostatevalue',arcidentifier,symkernel);
        optimalcostatename=fieldnames(optimalcostatevalue.value);
        if ~isempty(explicitcontroldynamicsindex.value)
            canonicalsystem=retrievemodelinformation(ocStruct,'canonicalsystem','',symkernel);
            invLSM=mycell2sym(invoptimalcontroldynamicsleftside.value,'matrix');
            RSM1=mycell2sym(pontryaginfunctionDuDX.value,'matrix');
            if ~isempty(pontryaginfunctionDlmmcDX.value)
                RSM2=mycell2sym(pontryaginfunctionDlmmcDX.value,'matrix');
            else
                RSM2=sym([]);
            end
            RSM=[RSM1;RSM2'];
            dxdt=mycell2sym(canonicalsystem.value,'vector');
            dudt=-invLSM*RSM*dxdt;
            if numel(explicitcontroldynamicsindex.value)==1
                conversiontype='vector';
            else
                conversiontype='matrix';
            end
            dudt=string2cell(char(dudt(explicitcontroldynamicsindex.value)),conversiontype);
            for ii=1:numel(optimalcostatename)
                for jj=1:numel(dudt)
                    %dudt{jj}=ocmatsimple(ocmatsubs(dudt{jj},[optimalcostatename{ii} '=' optimalcostatevalue.value.(optimalcostatename{ii})]),symkernel);
                    dudt{jj}=ocmatsubs(dudt{jj},[optimalcostatename{ii} '=' optimalcostatevalue.value.(optimalcostatename{ii})],symkernel);
                end
            end
            data.value=dudt;
        else
            data.value='';
        end
        data.type='mathcellchar';
        data.description='inverse of the left side matrix of the optimal control dynamics';

    case 'tensoroptimalcontroldynamicsleftside'
        controlname=retrievemodelinformation(ocStruct,'controlname');
        controlvector=cell2vectorstring(controlname.value);
        zerolmmcindex=retrievemodelinformation(ocStruct,'zerolmmcindex',arcidentifier);
        lagrangemultipliercontrolname=retrievemodelinformation(ocStruct,'lagrangemultipliercontrolname');
        lagrangemultipliercontrolname.value(zerolmmcindex.value)=[];
        lagrangemultipliervector=cell2vectorstring(lagrangemultipliercontrolname.value);
        pontryaginfunction=retrievemodelinformation(ocStruct,'pontryaginfunction');
        pontryaginfunction4arc=retrievemodelinformation(ocStruct,'pontryaginfunction4arc',arcidentifier,symkernel);
        pontryaginfunctionDu2=mystr2sym(ocmathessian(pontryaginfunction.value,controlvector,symkernel));
        if ~isempty(lagrangemultipliervector)
            pontryaginfunctionDlmmcDu=sym(ocmatjacobian(removematrixstring(ocmatjacobian(['[' pontryaginfunction4arc.value ']'],controlvector,symkernel)),lagrangemultipliervector,symkernel));
        else
            pontryaginfunctionDlmmcDu=sym([]);
        end
        if ~isempty(lagrangemultipliercontrolname.value)
            optimalcontroldynamicsleftside=char([pontryaginfunctionDu2 pontryaginfunctionDlmmcDu;pontryaginfunctionDlmmcDu.' sym(zeros(length(lagrangemultipliercontrolname.value)))]);
        else
            optimalcontroldynamicsleftside=char(pontryaginfunctionDu2);
        end
        optimalvalue=retrievemodelinformation(ocStruct,'optimalvalue',arcidentifier);
        maximizingvariable=retrievemodelinformation(ocStruct,'maximizingvariable');
        equationvariablename=retrievemodelinformation(ocStruct,'equationvariablenameimplicit',arcidentifier);
        for ii=1:numel(maximizingvariable.value)
            optimalcontroldynamicsleftside=ocmatsubs(optimalcontroldynamicsleftside,[maximizingvariable.value{ii} '=' optimalvalue.value.(maximizingvariable.value{ii})],symkernel);
        end
        for ii=1:numel(equationvariablename.value)
            T{ii}=string2cell(removematrixstring(ocmatdiff(optimalcontroldynamicsleftside,equationvariablename.value{ii},symkernel)),'matrix');
        end
        data.value=T;
        data.type='mathcellchar';
        data.description='left side matrix of the optimal control dynamics';

    case 'jacobianoptimalcontroldynamicsrightside'
        controlname=retrievemodelinformation(ocStruct,'controlname');
        controlvector=cell2vectorstring(controlname.value);
        statename=retrievemodelinformation(ocStruct,'statename');
        costatename=retrievemodelinformation(ocStruct,'costatename');
        statevector=cell2vectorstring([statename.value costatename.value]);
        specificstatedynamics=retrievemodelinformation(ocStruct,'specificstatedynamics',arcidentifier,'maple');
        specificadjointsystem=retrievemodelinformation(ocStruct,'specificadjointsystem',arcidentifier,'maple');
        dxdt=mycell2sym([specificstatedynamics.value specificadjointsystem.value],'vector');
        pontryaginfunction=retrievemodelinformation(ocStruct,'pontryaginfunction');
        J=removematrixstring(ocmatjacobian(['[' pontryaginfunction.value ']'],controlvector,symkernel));
        J=sym(ocmatjacobian(J,statevector,symkernel));
        RS=char(J*dxdt);
        optimalvalue=retrievemodelinformation(ocStruct,'optimalvalue',arcidentifier);
        maximizingvariable=retrievemodelinformation(ocStruct,'maximizingvariable');
        equationvariablename=retrievemodelinformation(ocStruct,'equationvariablenameimplicit',arcidentifier);
        vector=cell2vectorstring(equationvariablename.value);
        for ii=1:numel(maximizingvariable.value)
            RS=ocmatsubs(RS,[maximizingvariable.value{ii} '=' optimalvalue.value.(maximizingvariable.value{ii})],symkernel);
        end
        T=string2cell(char(jacobian(sym(RS),vector)),'charmatrix');
        data.value=T;
        data.type='mathcellchar';
        data.description='left side matrix of the optimal control dynamics';

    case {'pontryaginfunctionDuDX','pontryaginfunctionDXDu'}
        try
            data.value=ocStruct.pontryaginfunction.derivative.DuDX.term;
        catch
            controlname=retrievemodelinformation(ocStruct,'controlname');
            controlvector=cell2vectorstring(controlname.value);
            statename=retrievemodelinformation(ocStruct,'statename');
            costatename=retrievemodelinformation(ocStruct,'costatename');
            statevector=cell2vectorstring([statename.value costatename.value]);
            pontryaginfunction=retrievemodelinformation(ocStruct,'pontryaginfunction');
            J=removematrixstring(ocmatjacobian(['[' pontryaginfunction.value ']'],controlvector,symkernel));
            data.value=string2cell(ocmatjacobian(J,statevector,symkernel),'matrix');
        end
        data.type='mathcellchar';
        data.description='right side of the optimal control dynamics (-Pux*dotx-Pul*dotl)';

    case {'pontryaginfunctionDuDlmmc','pontryaginfunctionDlmmcDu'}
        try
            data.value=ocStruct.pontryaginfunction.derivative.DuDlmmc.term;
        catch
            controlname=retrievemodelinformation(ocStruct,'controlname');
            controlvector=cell2vectorstring(controlname.value);
            lagrangemultipliercontrolname=retrievemodelinformation(ocStruct,'lagrangemultipliercontrolname');
            lagrangemultipliervector=cell2vectorstring(lagrangemultipliercontrolname.value);
            pontryaginfunction=retrievemodelinformation(ocStruct,'pontryaginfunction');
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
        statename=retrievemodelinformation(ocStruct,'statename');
        costatename=retrievemodelinformation(ocStruct,'costatename');
        statevector=cell2vectorstring([statename.value costatename.value]);
        zerolmmcindex=retrievemodelinformation(ocStruct,'zerolmmcindex',arcidentifier);
        lagrangemultipliercontrolname=retrievemodelinformation(ocStruct,'lagrangemultipliercontrolname');
        lagrangemultipliercontrolname.value(zerolmmcindex.value)=[];
        lagrangemultipliervector=cell2vectorstring(lagrangemultipliercontrolname.value);
        pontryaginfunction=retrievemodelinformation(ocStruct,'pontryaginfunction');
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
            controlname=retrievemodelinformation(ocStruct,'controlname');
            controlvector=cell2vectorstring(controlname.value);
            lagrangemultipliercontrolname=retrievemodelinformation(ocStruct,'lagrangemultipliercontrolname');
            lagrangemultipliervector=cell2vectorstring(lagrangemultipliercontrolname.value);
            pontryaginfunction=retrievemodelinformation(ocStruct,'pontryaginfunction');
            pontryaginfunctionDuDlmmc=removematrixstring(ocmatjacobian(['[' pontryaginfunction.value ']'],controlvector,symkernel));
            pontryaginfunctionDuDlmmc=ocmatjacobian(pontryaginfunctionDuDlmmc,lagrangemultipliervector,symkernel);
        end
        zerolmmcindex=retrievemodelinformation(ocStruct,'zerolmmcindex',arcidentifier);
        pontryaginfunctionDuDlmmc(:,zerolmmcindex.value)=[];
        data.value=string2cell(removematrixstring(char(pontryaginfunctionDuDlmmc)),'charmatrix');
        data.type='mathcellchar';
        data.description='right side of the optimal control dynamics (-Pux*dotx-Pul*dotl)';

    case {'zeros4optimalcontroldynamics'}
        zerolmmcindex=retrievemodelinformation(ocStruct,'zerolmmcindex',arcidentifier);
        lagrangemultipliercontrolname=retrievemodelinformation(ocStruct,'lagrangemultipliercontrolname');
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
            maximizingimplicitvariable=retrievemodelinformation(ocStruct,'maximizingimplicitvariable',arcidentifier);
            controlname=retrievemodelinformation(ocStruct,'controlname');
            zerolmmcindex=retrievemodelinformation(ocStruct,'zerolmmcindex',arcidentifier);
            lagrangemultipliercontrolname=retrievemodelinformation(ocStruct,'lagrangemultipliercontrolname');
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
        pontryaginfunction=retrievemodelinformation(ocStruct,'pontryaginfunction');
        zerolmmcindex=retrievemodelinformation(ocStruct,'zerolmmcindex',arcidentifier);
        zerolmscindex=retrievemodelinformation(ocStruct,'zerolmscindex',arcidentifier);
        lagrangemultipliercontrolname=retrievemodelinformation(ocStruct,'lagrangemultipliercontrolname');
        lagrangemultiplierstatename=retrievemodelinformation(ocStruct,'lagrangemultiplierstatename');
        lmmcequalzero='';
        lmscequalzero='';
        ii=0;
        for idx=zerolmmcindex.value
            ii=ii+1;
            lmmcequalzero{ii}=[lagrangemultipliercontrolname.value{idx} '=0'];
        end
        zerolmmcvector=cell2vectorstring(lmmcequalzero);
        ii=0;
        for idx=zerolmscindex.value
            ii=ii+1;
            lmscequalzero{ii}=[lagrangemultiplierstatename.value{idx} '=0'];
        end
        zerolmscvector=cell2vectorstring(lmscequalzero);
        if ~isempty(zerolmmcvector)
            pontryaginfunction.value=ocmatsubs(pontryaginfunction.value,zerolmmcvector,symkernel);
        end
        if ~isempty(zerolmscvector)
            pontryaginfunction.value=ocmatsubs(pontryaginfunction.value,zerolmscvector,symkernel);
        end
        data.value=pontryaginfunction.value;
        data.type='mathchar';
        data.description='specific Hamilton function for actual constraint combination';

    case 'maximizingsolution'
        % Pontryagin's maximumprinciple is used to derive explicit formulas
        % of control and Lagrange multiplier values satisfying the first
        % order necessary conditions
        solution=[];
        totalmaximizingvariable=retrievemodelinformation(ocStruct,'maximizingvariable');
        maximizingexplicitvariable=retrievemodelinformation(ocStruct,'maximizingexplicitvariable',arcidentifier);
        maximizingimplicitvariable=retrievemodelinformation(ocStruct,'maximizingimplicitvariable',arcidentifier);
        zerolmmc=retrievemodelinformation(ocStruct,'zerolmmc',arcidentifier);
        nonzerolmmc=retrievemodelinformation(ocStruct,'nonzerolmmc',arcidentifier);
        zerolmsc=retrievemodelinformation(ocStruct,'zerolmsc',arcidentifier);
        maximizingderivativevariable=retrievemodelinformation(ocStruct,'maximizingderivativevariable',arcidentifier);
        controlname=retrievemodelinformation(ocStruct,'controlname');
        lagrangemultipliercontrolname=retrievemodelinformation(ocStruct,'lagrangemultipliercontrolname');
        lagrangemultiplierstatename=retrievemodelinformation(ocStruct,'lagrangemultiplierstatename');
        controlnum=retrievemodelinformation(ocStruct,'controlnum');
        inequalitycontrolconstraintnum=retrievemodelinformation(ocStruct,'inequalitycontrolconstraintnum');
        inequalitystateconstraintnum=retrievemodelinformation(ocStruct,'inequalitystateconstraintnum');
        transform2statecontrolspace=retrievemodelinformation(ocStruct,'transform2statecontrolspace',arcidentifier);
        pontryaginfunction4arc=retrievemodelinformation(ocStruct,'pontryaginfunction4arc',arcidentifier,symkernel);
        idx1=ocmatrank(pontryaginfunction4arc.value,cell2vectorstring(nonzerolmmc.value),cell2vectorstring(maximizingimplicitvariable.value));
        if ~isempty(idx1)
            nonzerolmmc.value(idx1)=[];
            maximizingvariable.value=[controlname.value nonzerolmmc.value];
        else
            maximizingvariable.value=maximizingexplicitvariable.value;
        end
        idx=ocmatrank(pontryaginfunction4arc.value,cell2vectorstring(maximizingvariable.value),cell2vectorstring(maximizingexplicitvariable.value));

        if transform2statecontrolspace.value
            controlname4controlstate=retrievemodelinformation(ocStruct,'controlname4controlstate',arcidentifier);
            optimalcostatevalue=retrievemodelinformation(ocStruct,'optimalcostatevalue',arcidentifier,symkernel);
            costatename=retrievemodelinformation(ocStruct,'costatename');
            explicitcontroldynamicsindex=retrievemodelinformation(ocStruct,'explicitcontroldynamicsindex',arcidentifier);
            maximizingexplicitvariable.value(explicitcontroldynamicsindex.value)=[];
            maximizingderivativevariable.value(explicitcontroldynamicsindex.value)=[];
        end
        if ~isempty(maximizingexplicitvariable.value)
            equationvariable=cell2vectorstring([maximizingexplicitvariable.value]);
            derivativevariable=cell2vectorstring([maximizingvariable.value(idx)]);
            equation=removematrixstring(ocmatjacobian(['[' pontryaginfunction4arc.value ']'],derivativevariable,symkernel));
            if transform2statecontrolspace.value
                for ii=1:length(costatename.value)
                    equation=removematrixstring(ocmatsubs(equation,[costatename.value{ii} '=' optimalcostatevalue.value.(costatename.value{ii})],symkernel));
                end
            end
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
                    for ii=1:numel(zerolmsc.value)
                        % the Lagrange multipliers being zero (inactive
                        % constraint) are added
                        solution(jj).(zerolmsc.value{ii})='0';
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
                for ii=1:numel(zerolmsc.value)
                    % the Lagrange multipliers being zero (inactive
                    % constraint) are added
                    solution(jj).(zerolmsc.value{ii})='0';
                end
            end
            % the order of the solution fields are controlname,
            % lmmcname
            if transform2statecontrolspace.value
                for ii=1:numel(explicitcontroldynamicsindex.value)
                    for jj=1:numel(solution)
                        solution(jj).(controlname4controlstate.value{ii})=controlname4controlstate.value{ii};
                    end
                end
            end
        elseif ~isempty(maximizingimplicitvariable.value)
            for ii=1:numel(maximizingimplicitvariable.value)
                solution.(maximizingimplicitvariable.value{ii})=maximizingimplicitvariable.value{ii};
            end
            for ii=1:numel(zerolmmc.value)
                % the Lagrange multipliers being zero (inactive
                % constraint) are added
                solution.(zerolmmc.value{ii})='0';
            end
            for ii=1:numel(zerolmsc.value)
                % the Lagrange multipliers being zero (inactive
                % constraint) are added
                solution.(zerolmsc.value{ii})='0';
            end
        else
            if transform2statecontrolspace.value
                for ii=1:numel(explicitcontroldynamicsindex.value)
                    for jj=1:max(1,numel(solution))
                        solution(jj).(controlname4controlstate.value{ii})=controlname4controlstate.value{ii};
                    end
                end
            else
                solution=[];
            end
        end
        solution=orderfields(solution,totalmaximizingvariable.value);
        for ii=1:numel(solution)
            for jj=1:controlnum.value
                solution(ii).control{jj}=solution(ii).(controlname.value{jj});
            end
            for jj=1:inequalitycontrolconstraintnum.value
                solution(ii).lagrangemultcc{jj}=solution(ii).(lagrangemultipliercontrolname.value{jj});
            end
            for jj=1:inequalitystateconstraintnum.value
                solution(ii).lagrangemultsc{jj}=solution(ii).(lagrangemultiplierstatename.value{jj});
            end
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
            pontryaginfunctionDx=retrievemodelinformation(ocStruct,'pontryaginfunctionDx','',symkernel);
            statenum=retrievemodelinformation(ocStruct,'statenum');
            costatename=retrievemodelinformation(ocStruct,'costatename');
            discountratevariable=retrievemodelinformation(ocStruct,'discountratevariable');
            for ii=1:statenum.value
                data.value{ii}=[discountratevariable.value '*' costatename.value{ii} '-(' pontryaginfunctionDx.value{ii} ')'];
            end
        end
        data.type='mathcellchar';
        data.description='costate dynamics';
    case 'adjointsystemcont'
        % the adjoint system, where each cell contains one costate dynamics
        arc=arcidentifier2field(arcidentifier);
        if isfield(ocStruct,'foc') && isfield(ocStruct.foc,'adjointsystem') ...
                && isfield(ocStruct.foc.adjointsystem,'dynamics') ...
                && isfield(ocStruct.foc.adjointsystem.dynamics,arc) ...
                && isfield(ocStruct.foc.adjointsystem.dynamics.(arc),'ode')
            data.value=ocStruct.foc.adjointsystem.dynamics.(arc).ode.term;
        else
            pontryaginfunctionDx=retrievemodelinformation(ocStruct,'pontryaginfunctionDx','',symkernel);
            statenum=retrievemodelinformation(ocStruct,'statenum');
            costatename=retrievemodelinformation(ocStruct,'costatename');
            discountratevariable=retrievemodelinformation(ocStruct,'discountratevariable');
            for ii=1:statenum.value
                data.value{ii}=[discountratevariable.value '*' costatename.value{ii} '-contpar*(' pontryaginfunctionDx.value{ii} ')'];
            end
        end
        data.type='mathcellchar';
        data.description='costate dynamics';

%     case 'algebraicequation'
%         % if the Hamiltonian maximizing condition does not allow an
%         % explicit algebraic expression for the control values, the
%         % corresponding equations are added as algebaric equations yielding
%         % a system of algebraic differential equations
%         arc=arcidentifier2field(arcidentifier);
%         %if false
%         if isfield(ocStruct,'foc') && isfield(ocStruct.foc,'adjointsystem') ...
%                 && isfield(ocStruct.foc.adjointsystem,'algebraicequation') ...
%                 && isfield(ocStruct.foc.adjointsystem.algebraicequation,arc)
%             data.value=ocStruct.foc.adjointsystem.algebraicequation.(arc).term;
%         else
%             pontryaginfunction4arc=retrievemodelinformation(ocStruct,'pontryaginfunction4arc',arcidentifier,symkernel);
%             maximizingimplicitderivativevariable=retrievemodelinformation(ocStruct,'maximizingimplicitderivativevariable',arcidentifier);
%             if ~isempty(maximizingimplicitderivativevariable.value)
%                 derivativevariable=cell2vectorstring([maximizingimplicitderivativevariable.value]);
%                 data.value=string2cell(removematrixstring(ocmatjacobian(['[' pontryaginfunction4arc.value ']'],derivativevariable,symkernel)),'vector');
%             end
%         end
%         data.type='mathcellchar';
%         data.description='equation to determine the implicitly given control values';

    case {'algebraicequation','algebraicequationimplicit'}
        % if the Hamiltonian maximizing condition does not allow an
        % explicit algebraic expression for the control values, the
        % corresponding equations are added as algebraic equations yielding
        % a system of algebraic differential equations
        %maximizingimplicitvariable=retrievemodelinformation(ocStruct,'maximizingimplicitvariable',arcidentifier);
        maximizingexplicitvariable=retrievemodelinformation(ocStruct,'maximizingexplicitvariable',arcidentifier);
        maximizingimplicitvariable=retrievemodelinformation(ocStruct,'maximizingimplicitvariable',arcidentifier);
        nonzerolmmc=retrievemodelinformation(ocStruct,'nonzerolmmc',arcidentifier);
        pontryaginfunction4arc=retrievemodelinformation(ocStruct,'pontryaginfunction4arc',arcidentifier,symkernel);
        totalvariable=maximizingexplicitvariable.value;
        idx1=ocmatrank(pontryaginfunction4arc.value,cell2vectorstring(nonzerolmmc.value),cell2vectorstring(maximizingimplicitvariable.value));
        if ~isempty(idx1)
            maximizingvariable=nonzerolmmc.value(idx1);
        else
            maximizingvariable=maximizingimplicitvariable.value;
        end
        %maximizingvariable=[controlname.value nonzerolmmc.value];

%         idx=ocmatrank(pontryaginfunction4arc.value,cell2vectorstring(maximizingvariable),cell2vectorstring(totalvariable));
%         totalvariable=maximizingvariable(setdiff(1:length(maximizingvariable),idx));
%         if ~isempty(nonzerolmmc.value)
%             remidx=[];
%             for ii=1:length(nonzerolmmc.value)
%                 remidx=[remidx find(strcmp(maximizingexplicitvariable.value,nonzerolmmc.value{ii}))];
%             end
%             maximizingexplicitvariable.value(remidx)=[];
%             derivativevariable=cell2vectorstring(nonzerolmmc.value);
%             mcu=ocmatjacobian(['[' pontryaginfunction4arc.value ']'],derivativevariable,symkernel);
%             derivativevariable=cell2vectorstring(maximizingexplicitvariable.value);
%             mcu=ocmatjacobian(mcu,derivativevariable,symkernel);
%             rk=double(rank(sym(mcu)));
%             maximizingexplicitvariable.value(1:rk)=[];
%         end
%         remidx=[];
%         for ii=1:length(maximizingimplicitderivativevariable.value)
%             remidx=[remidx find(strcmp(totalimplicitvariable.value,maximizingimplicitderivativevariable.value{ii}))];
%         end
%         totalimplicitvariable.value(remidx)=[];
        
        %derivativevariable=cell2vectorstring([maximizingimplicitderivativevariable.value maximizingexplicitvariable.value nonzerolmmc.value]);
        derivativevariable=cell2vectorstring(maximizingvariable);
        if ~isempty(derivativevariable)
            %derivativevariable=cell2vectorstring([totalimplicitvariable.value]);
            data.value=string2cell(removematrixstring(ocmatjacobian(['[' pontryaginfunction4arc.value ']'],derivativevariable,symkernel)),'vector');
        end
        data.type='mathcellchar';
        data.description='equation to determine the implicitly given control values';

    case 'specificstatedynamics'
        % the canonical system consisting of the state dynamics, adjoint
        % system and possible algebraic equations for implicitly given
        % control values.
        statedynamics=retrievemodelinformation(ocStruct,'statedynamics');
        optimalvalue=retrievemodelinformation(ocStruct,'optimalvalue',arcidentifier);
        maximizingvariable=retrievemodelinformation(ocStruct,'maximizingvariable');
        dxdt=statedynamics.value;

        for ii=1:numel(maximizingvariable.value)
            for jj=1:numel(dxdt)
                dxdt{jj}=ocmatsubs(dxdt{jj},[maximizingvariable.value{ii} '=' optimalvalue.value.(maximizingvariable.value{ii})],symkernel);
            end
        end
        data.value=dxdt;
        data.type='mathcellchar';
        data.description='state dynamics';

    case 'specificadjointsystem'
        % the canonical system consisting of the state dynamics, adjoint
        % system and possible algebraic equations for implicitly given
        % control values.
        adjointsystem=retrievemodelinformation(ocStruct,'adjointsystem','',symkernel);
        optimalvalue=retrievemodelinformation(ocStruct,'optimalvalue',arcidentifier);
        maximizingvariable=retrievemodelinformation(ocStruct,'maximizingvariable');
        dldt=adjointsystem.value;

        for ii=1:numel(maximizingvariable.value)
            for jj=1:numel(dldt)
                dldt{jj}=ocmatsubs(dldt{jj},[maximizingvariable.value{ii} '=' optimalvalue.value.(maximizingvariable.value{ii})],symkernel);
            end
        end
        data.value=dldt;
        data.type='mathcellchar';
        data.description='costate dynamics';

    case 'specificcanonicalsystem'
        % the canonical system consisting of the state dynamics, adjoint
        % system and possible algebraic equations for implicitly given
        % control values.
        transform2statecontrolspace=retrievemodelinformation(ocStruct,'transform2statecontrolspace',arcidentifier);
        if ~transform2statecontrolspace.value
            adjointsystem=retrievemodelinformation(ocStruct,'adjointsystem','',symkernel);
            algebraicequation=retrievemodelinformation(ocStruct,'algebraicequation',arcidentifier,symkernel);
            statedynamics=retrievemodelinformation(ocStruct,'statedynamics');
            optimalvalue=retrievemodelinformation(ocStruct,'optimalvalue',arcidentifier);
            maximizingvariable=retrievemodelinformation(ocStruct,'maximizingvariable');
            inequalitystateconstraintnum=retrievemodelinformation(ocStruct,'inequalitystateconstraintnum');
            dxdt=[statedynamics.value adjointsystem.value algebraicequation.value];
            for ii=1:numel(maximizingvariable.value)
                for jj=1:numel(dxdt)
                    if ~strcmp(dxdt{jj},'0')%~(dxdt{jj}==sym('0'))
                        dxdt{jj}=ocmatsubs(dxdt{jj},[maximizingvariable.value{ii} '=' optimalvalue.value.(maximizingvariable.value{ii})],symkernel);
                    end
                end
            end
            if inequalitystateconstraintnum.value
                explicitstatevalue=retrievemodelinformation(ocStruct,'explicitstatevalue',arcidentifier,symkernel);
                statename=retrievemodelinformation(ocStruct,'statename');
                for ii=1:length(statename.value)
                    for jj=1:numel(dxdt)
                        if ~strcmp(dxdt{jj},'0')%~(dxdt{jj}==sym('0'))
                            dxdt{jj}=ocmatsubs(dxdt{jj},[statename.value{ii} '=' explicitstatevalue.value{ii}],symkernel);
                        end
                    end
                end
            end
%             if inequalitystateconstraintnum.value
%                 nonzerolmscindex=retrievemodelinformation(ocStruct,'nonzerolmscindex',arcidentifier);
%                 if ~isempty(nonzerolmscindex.value)
%                     statename=retrievemodelinformation(ocStruct,'statename');
%                     explicitstatevalue=retrievemodelinformation(ocStruct,'explicitstatevalue',arcidentifier,symkernel);
%                     for ii=1:numel(statename.value)
%                         for jj=1:numel(dxdt)
%                             dxdt{jj}=ocmatsubs(dxdt{jj},[statename.value{ii} '=' explicitstatevalue.value{ii}]);
%                         end
%                     end
%                 end
%             end
            data.value=dxdt;
            data.type='mathcellchar';
            data.description='state dynamics, adjoint system (and algebraic equations for implicit control values)';
        else
            data=retrievemodelinformation(ocStruct,'specificcanonicalsysteminstatecontrol',arcidentifier,symkernel);
        end
        
    case 'equilibriumequation'
        specificcanonicalsystem=retrievemodelinformation(ocStruct,'specificcanonicalsystem',arcidentifier,symkernel);
        explicitstateindex=retrievemodelinformation(ocStruct,'explicitstateindex',arcidentifier);
        if ~isempty(explicitstateindex.value)
            inequalitystateconstraint=retrievemodelinformation(ocStruct,'inequalitystateconstraint');
            nonzerolmscindex=retrievemodelinformation(ocStruct,'nonzerolmscindex',arcidentifier);
            for ii=1:length(explicitstateindex.value)
                specificcanonicalsystem.value{explicitstateindex.value(ii)}=inequalitystateconstraint.value{nonzerolmscindex.value(ii)};
            end
        end
        data.value=specificcanonicalsystem.value;
        data.type='mathcellchar';
        data.description='';

    case 'optimalvalue4statecontrol'
        explicitcontroldynamicsindex=retrievemodelinformation(ocStruct,'explicitcontroldynamicsindex',arcidentifier);
        controlname=retrievemodelinformation(ocStruct,'controlname');
        explicitcontroldynamicsname=controlname.value(explicitcontroldynamicsindex.value);
        optimalvalue=retrievemodelinformation(ocStruct,'optimalvalue',arcidentifier);
        for ii=1:numel(explicitcontroldynamicsindex.value)
            optimalvalue.value.(explicitcontroldynamicsname{ii})=explicitcontroldynamicsname{ii};
        end
        data.value=optimalvalue.value;
        data.type='mathcellchar';
        data.description='state dynamics, adjoint system (and algebraic equations for implicit control values)';

    case 'specificcanonicalsysteminstatecontrol'
        % the canonical system consisting of the state dynamics, adjoint
        % system and possible algebraic equations for implicitly given
        % control values.
        explicitoptimalcontroldynamics=retrievemodelinformation(ocStruct,'explicitoptimalcontroldynamics',arcidentifier,symkernel);
        optimalcostatevalue=retrievemodelinformation(ocStruct,'optimalcostatevalue',arcidentifier,symkernel);
        replacedcostatename=retrievemodelinformation(ocStruct,'replacedcostatename',arcidentifier);
        optimalcostatename=replacedcostatename.value;
        replacedcostateindex=retrievemodelinformation(ocStruct,'replacedcostateindex',arcidentifier);
        adjointsystem=retrievemodelinformation(ocStruct,'adjointsystem','',symkernel);
        adjointsystem.value(replacedcostateindex.value)=explicitoptimalcontroldynamics.value;
        algebraicequation=retrievemodelinformation(ocStruct,'algebraicequation',arcidentifier,symkernel);
        inequalitystateconstraintnum=retrievemodelinformation(ocStruct,'inequalitystateconstraintnum');
        statedynamics=retrievemodelinformation(ocStruct,'statedynamics');
        optimalvalue4statecontrol=retrievemodelinformation(ocStruct,'optimalvalue4statecontrol',arcidentifier);
        maximizingvariable=retrievemodelinformation(ocStruct,'maximizingvariable');
        dxdt=[statedynamics.value adjointsystem.value algebraicequation.value];
        for ii=1:numel(optimalcostatename)
            for jj=1:numel(dxdt)
                %dudt{jj}=ocmatsimple(ocmatsubs(dudt{jj},[optimalcostatename{ii} '=' optimalcostatevalue.value.(optimalcostatename{ii})]),symkernel);
                dxdt{jj}=ocmatsubs(dxdt{jj},[optimalcostatename{ii} '=' optimalcostatevalue.value.(optimalcostatename{ii})],symkernel);
            end
        end
        for ii=1:numel(maximizingvariable.value)
            for jj=1:numel(dxdt)
                dxdt{jj}=ocmatsubs(dxdt{jj},[maximizingvariable.value{ii} '=' optimalvalue4statecontrol.value.(maximizingvariable.value{ii})],symkernel);
            end
        end
        if inequalitystateconstraintnum.value
            nonzerolmscindex=retrievemodelinformation(ocStruct,'nonzerolmscindex',arcidentifier);
            if ~isempty(nonzerolmscindex.value)
                statename=retrievemodelinformation(ocStruct,'statename');
                explicitstatevalue=retrievemodelinformation(ocStruct,'explicitstatevalue',arcidentifier,symkernel);
                for ii=1:numel(statename.value)
                    for jj=1:numel(dxdt)
                        dxdt{jj}=ocmatsubs(dxdt{jj},[statename.value{ii} '=' explicitstatevalue.value{ii}],symkernel);
                    end
                end
            end
        end
        data.value=dxdt;
        data.type='mathcellchar';
        data.description='state dynamics, adjoint system (and algebraic equations for implicit control values)';

    case 'specificpontryaginfunction'
        % the canonical system consisting of the state dynamics, adjoint
        % system and possible algebraic equations for implicitly given
        % control values.
        pontryaginfunction=retrievemodelinformation(ocStruct,'pontryaginfunction','',symkernel);
        optimalvalue=retrievemodelinformation(ocStruct,'optimalvalue',arcidentifier,symkernel);
        inequalitystateconstraintnum=retrievemodelinformation(ocStruct,'inequalitystateconstraintnum');
        maximizingvariable=retrievemodelinformation(ocStruct,'maximizingvariable');
        P=pontryaginfunction.value;
        for ii=1:numel(maximizingvariable.value)
            P=ocmatsubs(P,[maximizingvariable.value{ii} '=' optimalvalue.value.(maximizingvariable.value{ii})],symkernel);
        end
        if inequalitystateconstraintnum.value
            nonzerolmscindex=retrievemodelinformation(ocStruct,'nonzerolmscindex',arcidentifier);
            if ~isempty(nonzerolmscindex.value)
                statename=retrievemodelinformation(ocStruct,'statename');
                explicitstatevalue=retrievemodelinformation(ocStruct,'explicitstatevalue',arcidentifier,symkernel);
                for ii=1:numel(statename.value)
                    P=ocmatsubs(P,[statename.value{ii} '=' explicitstatevalue.value{ii}],symkernel);
                end
            end
        end
        data.value=P;
        data.type='mathcellchar';
        data.description='state dynamics, adjoint system (and algebraic equations for implicit control values)';

    case 'canonicalsystem'
        % the canonical system consisting of the state dynamics, adjoint
        % system and possible algebraic equations for implicitly given
        % control values.
        adjointsystem=retrievemodelinformation(ocStruct,'adjointsystem','',symkernel);
        statedynamics=retrievemodelinformation(ocStruct,'statedynamics');
        dxdt=[statedynamics.value adjointsystem.value];
        data.value=dxdt;
        data.type='mathcellchar';
        data.description='state dynamics, adjoint system';
        
    case 'canonicalsystemcont'
        % the canonical system consisting of the state dynamics, adjoint
        % system and possible algebraic equations for implicitly given
        % control values.
        adjointsystem=retrievemodelinformation(ocStruct,'adjointsystemcont','',symkernel);
        statedynamics=retrievemodelinformation(ocStruct,'statedynamics');
        dxdt=[statedynamics.value adjointsystem.value];
        data.value=dxdt;
        data.type='mathcellchar';
        data.description='state dynamics, adjoint system';

    case 'generalcanonicalsystem'
        % the canonical system consisting of the state dynamics, adjoint
        % system and possible algebraic equations for implicitly given
        % control values.
        adjointsystem=retrievemodelinformation(ocStruct,'adjointsystem','',symkernel);
        algebraicequation=retrievemodelinformation(ocStruct,'algebraicequation',arcidentifier,symkernel);
        statedynamics=retrievemodelinformation(ocStruct,'statedynamics');
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
        canonicalsystem=retrievemodelinformation(ocStruct,'specificcanonicalsysteminstatecontrol',arcidentifier,symkernel);
        dxdt=cell2vectorstring(canonicalsystem.value);
        equationvariablename=retrievemodelinformation(ocStruct,'equationvariablename',arcidentifier);
        variable=cell2vectorstring(equationvariablename.value);
        data.value=string2cell(removematrixstring(ocmatjacobian(dxdt,variable,symkernel)),'matrix');
        data.type='mathcellchar';
        data.description='Jacobian of the canonical system with respect to state, costate and control/Lagrange multiplier variables, with explicit control dynamics';

    case 'statecostatejacobian'
        % Jacobian of the canonical system with respect to state, costate
        % and implicit control/Lagrange multiplier variables
        specificstatedynamics=retrievemodelinformation(ocStruct,'specificstatedynamics',arcidentifier,symkernel);
        specificadjointsystem=retrievemodelinformation(ocStruct,'specificadjointsystem',arcidentifier,symkernel);
        dxdt=cell2vectorstring([specificstatedynamics.value specificadjointsystem.value]);
        equationvariablename=retrievemodelinformation(ocStruct,'equationvariablename',arcidentifier);
        variable=cell2vectorstring(equationvariablename.value);
        data.value=string2cell(removematrixstring(ocmatjacobian(dxdt,variable,symkernel)),'matrix');
        data.type='mathcellchar';
        data.description='Jacobian of the canonical system with respect to state, costate and implicit control/Lagrange multiplier variables';

    case 'canonicalsystemhessian'
        % Jacobian of the canonical system with respect to state, costate
        % and implicit control/Lagrange multiplier variables
        canonicalsystem=retrievemodelinformation(ocStruct,'specificcanonicalsystem',arcidentifier,symkernel);
        equationvariablename=retrievemodelinformation(ocStruct,'equationvariablename',arcidentifier);
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
        canonicalsystem=retrievemodelinformation(ocStruct,'specificcanonicalsysteminstatecontrol',arcidentifier,symkernel);
        dxdt=cell2vectorstring(canonicalsystem.value);
        parametername=retrievemodelinformation(ocStruct,'parametername',arcidentifier);
        variable=cell2vectorstring(parametername.value);
        data.value=string2cell(removematrixstring(ocmatjacobian(dxdt,variable,symkernel)),'matrix');
        data.type='mathcellchar';
        data.description='Jacobian of the canonical system with respect to exogenous parameter variables';

    case 'canonicalsystemderivativetime'
        % Jacobian of the canonical system with respect to exogenous
        % parameter variables
        canonicalsystem=retrievemodelinformation(ocStruct,'specificcanonicalsystem',arcidentifier,symkernel);
        dxdt=cell2vectorstring(canonicalsystem.value);
        timevariable=retrievemodelinformation(ocStruct,'independent');
        variable=cell2vectorstring(timevariable.value);
        data.value=string2cell(removematrixstring(ocmatjacobian(dxdt,variable,symkernel)),'matrix');
        data.type='mathcellchar';
        data.description='Derivative of the canonical system with respect to the independent variable';

    case 'canonicalsystemtotalhessian'
        % Jacobian of the canonical system with respect to exogenous
        % parameter variables
        canonicalsystem=retrievemodelinformation(ocStruct,'specificcanonicalsystem',arcidentifier,symkernel);
        equationvariablename=retrievemodelinformation(ocStruct,'equationvariablename',arcidentifier);
        parametername=retrievemodelinformation(ocStruct,'parametername',arcidentifier);
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
        controlname=retrievemodelinformation(ocStruct,'controlname');
        controlnum=retrievemodelinformation(ocStruct,'controlnum');
        lagrangemultipliercontrolname=retrievemodelinformation(ocStruct,'lagrangemultipliercontrolname');
        inequalitycontrolconstraintnum=retrievemodelinformation(ocStruct,'inequalitycontrolconstraintnum');
        lagrangemultiplierstatename=retrievemodelinformation(ocStruct,'lagrangemultiplierstatename');
        inequalitystateconstraintnum=retrievemodelinformation(ocStruct,'inequalitystateconstraintnum');

        for ii=1:controlnum.value
            value.(controlname.value{ii})=ocStruct.foc.value.control.(arcidentifier2field(arcidentifier)).term{ii};
        end
        for ii=1:inequalitycontrolconstraintnum.value
            value.(lagrangemultipliercontrolname.value{ii})=ocStruct.foc.value.lagrangemultcc.(arcidentifier2field(arcidentifier)).term{ii};
        end
        for ii=1:inequalitystateconstraintnum.value
            value.(lagrangemultiplierstatename.value{ii})=ocStruct.foc.value.lagrangemultsc.(arcidentifier2field(arcidentifier)).term{ii};
        end
        data.type='structmathchar';
        data.value=value;
        data.description='optimal control values as derived from the Hamiltonian maximizing condition';

    case 'optimalcostatevalue'
        % the solutions from the Hamiltonian maximizing condition have
        % already to be stored in the ocStruct structure.
        explicitcontroldynamicsindex=retrievemodelinformation(ocStruct,'explicitcontroldynamicsindex',arcidentifier);
        replacedcostateindex=retrievemodelinformation(ocStruct,'replacedcostateindex',arcidentifier);
        controlname=retrievemodelinformation(ocStruct,'controlname');
        costatename=retrievemodelinformation(ocStruct,'costatename');
        costatenum=retrievemodelinformation(ocStruct,'costatenum');
        pontryaginfunction4arc=retrievemodelinformation(ocStruct,'pontryaginfunction4arc',arcidentifier,symkernel);
        replacedcostate=costatename.value(replacedcostateindex.value);
        equationvariable=cell2vectorstring([replacedcostate]);
        derivativevariable=cell2vectorstring([controlname.value(explicitcontroldynamicsindex.value)]);
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
            for ii=1:costatenum.value
                if ~isfield(solution,costatename.value{ii})
                    solution.(costatename.value{ii})=costatename.value{ii};
                end
            end
            solution=orderfields(solution,costatename.value);
        end
        data.type='structmathchar';
        data.value=solution;
        data.description='';



    case 'hamiltonianDu'
        % derivative of the Pontryaginfunction with respect to the control
        % variable(s)
        controlname=retrievemodelinformation(ocStruct,'controlname');
        controlvector=cell2vectorstring(controlname.value);
        hamiltonianfunction=retrievemodelinformation(ocStruct,'hamiltonianfunction');
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
        statename=retrievemodelinformation(ocStruct,'statename');
        statevector=cell2vectorstring(statename.value);
        salvagevalue=retrievemodelinformation(ocStruct,'salvagevalue');
        if isempty(salvagevalue.value)
            salvagevalue.value='0';
        end
        data.value=string2cell(ocmatjacobian(['[' salvagevalue.value ']'],statevector,symkernel),'vector');
        %        end
        data.type='mathcellchar';
        data.description='derivative of the Pontryaginfunction with respect to the state';
    otherwise
        ocmatmsg('No value returned. Property class ''%s'' is unknown.\n',propertyclass)
end



