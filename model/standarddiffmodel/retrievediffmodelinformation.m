function data=retrievediffmodelinformation(ocStruct,propertyclass,arcidentifier,symkernel)
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
        constraintcombination=retrievediffmodelinformation(ocStruct,'constraintcombination',arcidentifier);
        inequalitycontrolconstraintidentifier=retrievediffmodelinformation(ocStruct,'inequalitycontrolconstraintidentifier');
        inequalitycontrolconstraintnum=retrievediffmodelinformation(ocStruct,'inequalitycontrolconstraintnum');
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

    case 'statenamet'
        % variable name of the state(s) as defined by the user in the
        % initialization file
        independent=retrievediffmodelinformation(ocStruct,'independent');
        statename=retrievediffmodelinformation(ocStruct,'statename');
        data.type='cellchar';
        for ii=1:length(statename.value)
            data.value{ii}=[statename.value{ii} '_' independent.value];
        end
        data.description='name of state';

    case 'statenametp1'
        % variable name of the state(s) as defined by the user in the
        % initialization file
        independent=retrievediffmodelinformation(ocStruct,'independent');
        statename=retrievediffmodelinformation(ocStruct,'statename');
        data.type='cellchar';
        for ii=1:length(statename.value)
            data.value{ii}=[statename.value{ii} '_' independent.value 'p1'];
        end
        data.description='name of state';

    case 'costatename'
        % variable name of the costate(s) as defined by the user in the
        % initialization file
        data.type='cellchar';
        data.value=ocStruct.variable.costate.name;
        data.description='name of costate';

    case 'costatenamet'
        % variable name of the state(s) as defined by the user in the
        % initialization file
        independent=retrievediffmodelinformation(ocStruct,'independent');
        costatename=retrievediffmodelinformation(ocStruct,'costatename');
        data.type='cellchar';
        for ii=1:length(costatename.value)
            data.value{ii}=[costatename.value{ii} '_' independent.value];
        end
        data.description='name of costate';

    case 'costatenametp1'
        % variable name of the state(s) as defined by the user in the
        % initialization file
        independent=retrievediffmodelinformation(ocStruct,'independent');
        costatename=retrievediffmodelinformation(ocStruct,'costatename');
        data.type='cellchar';
        for ii=1:length(costatename.value)
            data.value{ii}=[costatename.value{ii} '_' independent.value 'p1'];
        end
        data.description='name of state';

    case 'controlname'
        % variable name of the control(s) as defined by the user in the
        % initialization file
        data.type='cellchar';
        data.value=ocStruct.variable.control.name;
        data.description='name of control';

    case 'controlnamet'
        % variable name of the state(s) as defined by the user in the
        % initialization file
        independent=retrievediffmodelinformation(ocStruct,'independent');
        controlname=retrievediffmodelinformation(ocStruct,'controlname');
        data.type='cellchar';
        for ii=1:length(controlname.value)
            data.value{ii}=[controlname.value{ii} '_' independent.value];
        end
        data.description='name of costate';

    case 'controlnametp1'
        % variable name of the state(s) as defined by the user in the
        % initialization file
        independent=retrievediffmodelinformation(ocStruct,'independent');
        controlname=retrievediffmodelinformation(ocStruct,'controlname');
        data.type='cellchar';
        for ii=1:length(controlname.value)
            data.value{ii}=[controlname.value{ii} '_' independent.value 'p1'];
        end
        data.description='name of state';

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

    case 'exogenousfunctionterm'
        % variable name of the control(s) as defined by the user in the
        % initialization file
        data.type='cellchar';
        exogenousfunctionname=retrievediffmodelinformation(ocStruct,'exogenousfunctionname');
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
        parametername=retrievediffmodelinformation(ocStruct,'parametername');
        data.value=zeros(1,numel(parametername.value));
        for ii=1:numel(parametername.value)
            data.value(ii)=ocStruct.parameter.variable.(parametername.value{ii});
        end
        data.description='vector of user provided parameter values';

    case 'discountrate'
        data.type='double';
        parametername=retrievediffmodelinformation(ocStruct,'parametername');
        parametervalue=retrievediffmodelinformation(ocStruct,'parametervalue');
        discountratevariable=retrievediffmodelinformation(ocStruct,'discountratevariable');
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

    case 'lagrangemultipliercontrolnamet'
        % variable name of the Lagrange multiplier(s) for inequality
        % (control) constraints as defined by the user in the
        % initialization file
        independent=retrievediffmodelinformation(ocStruct,'independent');
        lagrangemultipliercontrolname=retrievediffmodelinformation(ocStruct,'lagrangemultipliercontrolname');
        data.type='cellchar';
        for ii=1:length(lagrangemultipliercontrolname.value)
            data.value{ii}=[lagrangemultipliercontrolname.value{ii} '_' independent.value];
        end

    case 'lagrangemultipliercontrolnametp1'
        % variable name of the Lagrange multiplier(s) for inequality
        % (control) constraints as defined by the user in the
        % initialization file
        independent=retrievediffmodelinformation(ocStruct,'independent');
        lagrangemultipliercontrolname=retrievediffmodelinformation(ocStruct,'lagrangemultipliercontrolname');
        data.type='cellchar';
        for ii=1:length(lagrangemultipliercontrolname.value)
            data.value{ii}=[lagrangemultipliercontrolname.value{ii} '_' independent.value 'p1'];
        end

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
        maximizingimplicitderivativevariable=retrievediffmodelinformation(ocStruct,'maximizingimplicitderivativevariable',arcidentifier);
        data.type='integer';
        data.value=numel(maximizingimplicitderivativevariable.value);
        data.description='number of algebraic equations';

    case 'canonicalsystemequationnum'
        statenum=retrievediffmodelinformation(ocStruct,'statenum');
        costatenum=retrievediffmodelinformation(ocStruct,'costatenum');
        algebraicequationnum=retrievediffmodelinformation(ocStruct,'algebraicequationnum',arcidentifier);
        data.type='integer';
        data.value=statenum.value+costatenum.value+algebraicequationnum.value;
        data.description='total number of algebraic differential equations';

    case 'statedynamics'
        % ODEs describing the evolution of the state(s)
        data.type='mathchar';
        data.value=ocStruct.constraint.difference.state.term;
        data.description='ODEs of the state dynamics';

    case 'statedynamicsimplicit'
        % ODEs describing the evolution of the state(s)
        statedynamics=retrievediffmodelinformation(ocStruct,'statedynamics');
        statenametp1=retrievediffmodelinformation(ocStruct,'statenametp1');
        data.type='mathchar';
        for ii=1:numel(statedynamics.value)
            data.value{ii}=[statedynamics.value{ii} '-' statenametp1.value{ii}];
        end;
        data.description='ODEs of the state dynamics';

    case 'maptype'
        % ODEs describing the evolution of the state(s)
        data.type='mathchar';
        data.value=ocStruct.constraint.difference.maptype;
        data.description='ODEs of the state dynamics';

    case 'objectivetype'
        % type of the terms which build up the objective value function
        % (usually sum and Salvage value)
        data.type='mathchar';
        data.value=fieldnames(ocStruct.objective);
        data.description='type of terms which build up the objective value function';

    case 'objectivesummand'
        % summand for objective value function without the possible
        % discounting factor
        data.type='mathchar';
        data.value=ocStruct.objective.sum.function.term;
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
        if isfield(ocStruct.objective,'sum')
            data.type='char';
            data.value=ocStruct.objective.sum.discountrate;
            data.description='discount variable of the objective value function';
        end

    case 'discountfactor'
        % discount variable used in the discount factor of the objective
        % value function
        if isfield(ocStruct.objective,'sum')
            data.type='mathchar';
            data.value=ocStruct.objective.sum.discountfactor.term;
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
        objectivefunction=retrievediffmodelinformation(ocStruct,'objectivesummand');
        discountratevariable=retrievediffmodelinformation(ocStruct,'discountratevariable');
        independentvariable=retrievediffmodelinformation(ocStruct,'independent');
        if isfield(ocStruct.objective,'sum')
            data.type='mathchar';
            data.value=['exp(-' discountratevariable.value '*' independentvariable.value ')*(' objectivefunction.value ')'];
            data.description='the discounted objective summand';
        end

    case 'specificdiscobjectivefunction'
        discobjectivefunction=retrievediffmodelinformation(ocStruct,'discobjectivefunction');
        optimalvalue=retrievediffmodelinformation(ocStruct,'optimalvalue',arcidentifier);
        maximizingvariable=retrievediffmodelinformation(ocStruct,'maximizingvariable');
        for ii=1:numel(maximizingvariable.value)
            discobjectivefunction.value=ocmatsubs(discobjectivefunction.value,[maximizingvariable.value{ii} '=' optimalvalue.value.(maximizingvariable.value{ii})]);
        end
        data.type='mathchar';
        data.value=discobjectivefunction.value;
        data.description='the discounted objective summand';

    case 'discobjectivefunctionDX'
        %costatenamet=retrievediffmodelinformation(ocStruct,'costatenamet');
        statenamet=retrievediffmodelinformation(ocStruct,'statenamet');
        costatenametp1=retrievediffmodelinformation(ocStruct,'costatenametp1');
        %statenametp1=retrievediffmodelinformation(ocStruct,'statenametp1');
        statevector=cell2vectorstring([statenamet.value costatenametp1.value]);
        specificdiscobjectivefunction=retrievediffmodelinformation(ocStruct,'specificdiscobjectivefunction',arcidentifier);
        data.type='mathcellchar';
        data.value=string2cell(ocmatjacobian(['[' specificdiscobjectivefunction.value ']'],statevector,symkernel),'vector');
        data.description='derivative of the discounted objective summand with respect to the state and costate';

    case 'discobjectivefunctionDP'
        specificdiscobjectivefunction=retrievediffmodelinformation(ocStruct,'specificdiscobjectivefunction',arcidentifier);
        parametername=retrievediffmodelinformation(ocStruct,'parametername',arcidentifier);
        variable=cell2vectorstring(parametername.value);
        data.type='mathcellchar';
        data.value=string2cell(ocmatjacobian(['[' specificdiscobjectivefunction.value ']'],variable,symkernel),'vector');
        data.description='derivative of the discounted objective summand with respect to the state and costate';

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
        statename=retrievediffmodelinformation(ocStruct,'statename');
        costatename=retrievediffmodelinformation(ocStruct,'costatename');
        implicitcontrolname=retrievediffmodelinformation(ocStruct,'implicitnonlinearcontrol',arcidentifier);
        data.value=[statename.value costatename.value implicitcontrolname.value];
        data.description='dependent variables of the canonical system';

    case 'equationvariablenamet'
        % dependent variables of the canonical system for a specific arc
        data.type='cellchar';
        statenamet=retrievediffmodelinformation(ocStruct,'statenamet');
        costatenamet=retrievediffmodelinformation(ocStruct,'costatenamet');
        implicitnonlinearcontroltp1=retrievediffmodelinformation(ocStruct,'implicitnonlinearcontroltp1',arcidentifier);
        data.value=[statenamet.value costatenamet.value implicitnonlinearcontroltp1.value];
        data.description='dependent variables of the canonical system';

    case 'equationvariablenametp1'
        % dependent variables of the canonical system for a specific arc
        data.type='cellchar';
        statenametp1=retrievediffmodelinformation(ocStruct,'statenametp1');
        costatenametp1=retrievediffmodelinformation(ocStruct,'costatenametp1');
        implicitcontrolnametp1=retrievediffmodelinformation(ocStruct,'implicitnonlinearcontroltp1',arcidentifier);
        data.value=[statenametp1.value costatenametp1.value implicitcontrolnametp1.value];
        data.description='dependent variables of the canonical system';

    case 'explicitnonlinearcontroltp1'
        % explicit nonlinear control variables for a specific arc
        property=ocStruct.variable.control.(arcidentifier2field(arcidentifier)).property;
        independent=retrievediffmodelinformation(ocStruct,'independent');
        data.value=ocStruct.variable.control.name([property.explicit{:}]==1 & [property.linear{:}]==0);
        for ii=1:length(data.value)
            data.value{ii}=[data.value{ii} '_' independent.value 'p1'];
        end

        data.type='cellchar';
        data.description='explicit nonlinear control variable name';

    case 'implicitnonlinearcontrol'
        % implicit nonlinear control variables for a specific arc
        property=ocStruct.variable.control.(arcidentifier2field(arcidentifier)).property;
        data.value=ocStruct.variable.control.name([property.explicit{:}]==0 & [property.linear{:}]==0);
        data.type='cellchar';
        data.description='implicit nonlinear control variable name';

    case 'implicitnonlinearcontroltp1'
        % implicit nonlinear control variables for a specific arc
        property=ocStruct.variable.control.(arcidentifier2field(arcidentifier)).property;
        independent=retrievediffmodelinformation(ocStruct,'independent');
        data.value=ocStruct.variable.control.name([property.explicit{:}]==0 & [property.linear{:}]==0);
        for ii=1:length(data.value)
            data.value{ii}=[data.value{ii} '_' independent.value 'p1'];
        end
        data.type='cellchar';
        data.description='implicit nonlinear control variable name';

    case 'implicitnonlinearcontrolindex'
        % implicit nonlinear control variables for a specific arc
        property=ocStruct.variable.control.(arcidentifier2field(arcidentifier)).property;
        data.value=find([property.explicit{:}]==0 & [property.linear{:}]==0);
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
        controlnametp1=retrievediffmodelinformation(ocStruct,'controlnametp1');
        lagrangemultipliercontrolnametp1=retrievediffmodelinformation(ocStruct,'lagrangemultipliercontrolnametp1');
        data.value=[controlnametp1.value lagrangemultipliercontrolnametp1.value];
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
            explicitnonlinearcontroltp1=retrievediffmodelinformation(ocStruct,'explicitnonlinearcontroltp1',arcidentifier);
            nonzerolmmctp1=retrievediffmodelinformation(ocStruct,'nonzerolmmctp1',arcidentifier);
            data.value=[explicitnonlinearcontroltp1.value nonzerolmmctp1.value];
        end
        data.type='cellchar';
        data.description='variable names for the first order necessary conditions of the Hamiltonian maximizing condition';

    case 'maximizingimplicitderivativevariable'
        % variable names which are used to derive the first order necessary
        % conditions of the Hamiltonian maximizing condition and are only
        % calculated implicitly
        data=retrievediffmodelinformation(ocStruct,'maximizingvariable');
        maximizingderivativevariable=retrievediffmodelinformation(ocStruct,'maximizingderivativevariable',arcidentifier);
        zerolmmctp1=retrievediffmodelinformation(ocStruct,'zerolmmctp1',arcidentifier);
        totalexplicitvariable=[maximizingderivativevariable.value zerolmmctp1.value];
        for ii=1:numel(totalexplicitvariable)
            data.value(strcmp(data.value,totalexplicitvariable{ii}))=[];
        end
        data.description='implicit variable names for the first order necessary conditions of the Hamiltonian maximizing condition';

    case 'maximizingexplicitvariable'
        % variable names which are used to derive the first order necessary
        % conditions of the Hamiltonian maximizing condition and are
        % calculated explicitly
        arc=arcidentifier2field(arcidentifier);
        if isfield(ocStruct.foc,'generalinformation') && ...
                isfield(ocStruct.foc.generalinformation,arc) && ...
                isfield(ocStruct.foc.generalinformation.(arc),'maximizingexplicitvariable')
            data.value=ocStruct.foc.generalinformation.(arc).maximizingexplicitvariable.name;
        else
            explicitnonlinearcontroltp1=retrievediffmodelinformation(ocStruct,'explicitnonlinearcontroltp1',arcidentifier);
            nonzerolmmctp1=retrievediffmodelinformation(ocStruct,'nonzerolmmctp1',arcidentifier);
            data.value=[explicitnonlinearcontroltp1.value nonzerolmmctp1.value];
        end
        data.type='cellchar';
        data.description='explicit variable names for the first order necessary conditions of the Hamiltonian maximizing condition';

    case 'maximizingimplicitvariable'
        % variable names which are implicitly calculated by the Hamiltonian
        % maximizing condition
        data=retrievediffmodelinformation(ocStruct,'maximizingvariable');
        maximizingexplicitvariable=retrievediffmodelinformation(ocStruct,'maximizingexplicitvariable',arcidentifier);
        zerolmmctp1=retrievediffmodelinformation(ocStruct,'zerolmmctp1',arcidentifier);
        totalexplicitvariable=[maximizingexplicitvariable.value zerolmmctp1.value];
        for ii=1:numel(totalexplicitvariable)
            data.value(strcmp(data.value,totalexplicitvariable{ii}))=[];
        end
        data.description='implicit variable names of the Hamiltonian maximizing condition';
        
    case 'totalimplicitvariable'
        controlname=retrievediffmodelinformation(ocStruct,'controlname');
        totalimplicitnonlinearcontrolindex=retrievediffmodelinformation(ocStruct,'totalimplicitvariableindex');
        data.value=controlname.value(totalimplicitnonlinearcontrolindex.value);
        data.type='integer';
        data.description='index of the Lagrange multipliers corresponding to an inactive constraint';
        
    case 'totalimplicitvariableindex'
        arcidentifier=retrievediffmodelinformation(ocStruct,'arcidentifier');
        totalimplicitnonlinearcontrolindex=[];
        for ii=1:numel(arcidentifier.value)
            implicitnonlinearcontrolindex=retrievediffmodelinformation(ocStruct,'implicitnonlinearcontrolindex',arcidentifier.value{ii});
            totalimplicitnonlinearcontrolindex=[totalimplicitnonlinearcontrolindex implicitnonlinearcontrolindex.value];
        end
        data.value=unique(totalimplicitnonlinearcontrolindex);
        data.type='integer';
        data.description='index of the Lagrange multipliers corresponding to an inactive constraint';

    case 'equationvariablenameimplicit'
        % dependent variables of the canonical system for a specific arc
        data.type='cellchar';
        statename=retrievediffmodelinformation(ocStruct,'statename');
        costatename=retrievediffmodelinformation(ocStruct,'costatename');
        implicitcontrolname=retrievediffmodelinformation(ocStruct,'maximizingimplicitvariable',arcidentifier);
        data.value=[statename.value costatename.value implicitcontrolname.value];
        data.description='dependent variables of the canonical system';

    case 'zerolmmcindex'
        % returns the index of the Lagrange multipliers corresponding to an
        % inactive constraint for a specific arc
        inequalitycontrolconstraintidentifier=retrievediffmodelinformation(ocStruct,'inequalitycontrolconstraintidentifier');
        controlconstraintnum=retrievediffmodelinformation(ocStruct,'inequalitycontrolconstraintnum');
        constraintcombination=retrievediffmodelinformation(ocStruct,'constraintcombination',arcidentifier);
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
        lagrangemultipliercontrolname=retrievediffmodelinformation(ocStruct,'lagrangemultipliercontrolname');
        zerolmmcindex=retrievediffmodelinformation(ocStruct,'zerolmmcindex',arcidentifier);
        data.value=lagrangemultipliercontrolname.value(zerolmmcindex.value);
        data.type='cellchar';
        data.description='name of the Lagrange multipliers corresponding to an inactive constraint';

    case 'zerolmmctp1'
        % returns the index of the Lagrange multipliers corresponding to an
        % inactive constraint for a specific arc
        lagrangemultipliercontrolnametp1=retrievediffmodelinformation(ocStruct,'lagrangemultipliercontrolnametp1');
        zerolmmcindex=retrievediffmodelinformation(ocStruct,'zerolmmcindex',arcidentifier);
        data.value=lagrangemultipliercontrolnametp1.value(zerolmmcindex.value);
        data.type='cellchar';
        data.description='name of the Lagrange multipliers corresponding to an inactive constraint';

    case 'nonzerolmmc'
        % returns the name of the Lagrange multipliers corresponding to an
        % active constraint for a specific arc
        lagrangemultipliercontrolname=retrievediffmodelinformation(ocStruct,'lagrangemultipliercontrolname');
        nonzerolmmcindex=retrievediffmodelinformation(ocStruct,'nonzerolmmcindex',arcidentifier);
        data.value=lagrangemultipliercontrolname.value(nonzerolmmcindex.value);
        data.type='cellchar';
        data.description='name of the Lagrange multipliers corresponding to an active constraint';

    case 'nonzerolmmctp1'
        % returns the name of the Lagrange multipliers corresponding to an
        % active constraint for a specific arc
        lagrangemultipliercontrolnametp1=retrievediffmodelinformation(ocStruct,'lagrangemultipliercontrolnametp1');
        nonzerolmmcindex=retrievediffmodelinformation(ocStruct,'nonzerolmmcindex',arcidentifier);
        data.value=lagrangemultipliercontrolnametp1.value(nonzerolmmcindex.value);
        data.type='cellchar';
        data.description='name of the Lagrange multipliers corresponding to an active constraint';

    case 'nonzerolmmcindex'
        % returns the index of the Lagrange multipliers corresponding to an
        % active constraint for a specific arc
        inequalitycontrolconstraintnum=retrievediffmodelinformation(ocStruct,'inequalitycontrolconstraintnum');
        zerolmmcindex=retrievediffmodelinformation(ocStruct,'zerolmmcindex',arcidentifier);
        data.value=setdiff(1:inequalitycontrolconstraintnum.value,zerolmmcindex.value);
        data.type='integer';
        data.description='index of the Lagrange multipliers corresponding to an active constraint';

    case 'pontryaginfunction'
        % general Hamilton function including the "Lagrange terms" from
        % inequality contraints
        if isfield(ocStruct,'pontryaginfunction') && isfield(ocStruct.pontryaginfunction,'term')
            data.value=ocStruct.pontryaginfunction.term;
        else
            statenum=retrievediffmodelinformation(ocStruct,'statenum');
            costatenametp1=retrievediffmodelinformation(ocStruct,'costatenametp1');
            inequalitycontrolconstraint=retrievediffmodelinformation(ocStruct,'inequalitycontrolconstraint');
            inequalitycontrolconstraintnum=retrievediffmodelinformation(ocStruct,'inequalitycontrolconstraintnum');
            lagrangemultipliercontrolnametp1=retrievediffmodelinformation(ocStruct,'lagrangemultipliercontrolnametp1');
            
            % summand of the objective value function without discount factor
            g=retrievediffmodelinformation(ocStruct,'objectivesummand');
            % state dynamics
            xp1=retrievediffmodelinformation(ocStruct,'statedynamics');
            data.value=g.value;
            for ii=1:statenum.value
                data.value=[data.value '+' costatenametp1.value{ii} '*(' xp1.value{ii} ')'];
            end
            for ii=1:inequalitycontrolconstraintnum.value
                data.value=[data.value '+' lagrangemultipliercontrolnametp1.value{ii} '*(' inequalitycontrolconstraint.value{ii} ')'];
            end
        end
        data.type='mathchar';
        data.description='general Hamilton function P=H+mu*mc';

    case 'pontryaginfunctionDxt'
        % derivative of the Pontryaginfunction with respect to the state
        % variable(s)
        if isfield(ocStruct,'pontryaginfunction') && isfield(ocStruct.pontryaginfunction,'derivative') ...
                && isfield(ocStruct.pontryaginfunction.derivative,'Dxt')
            data.value=ocStruct.pontryaginfunction.derivative.Dxt.term;
        else
            statenamet=retrievediffmodelinformation(ocStruct,'statenamet');
            statevector=cell2vectorstring(statenamet.value);
            pontryaginfunction=retrievediffmodelinformation(ocStruct,'pontryaginfunction');
            data.value=string2cell(ocmatjacobian(['[' pontryaginfunction.value ']'],statevector,symkernel),'vector');
        end
        data.type='mathcellchar';
        data.description='derivative of the Pontryaginfunction with respect to the state';

    case 'pontryaginfunctionDX'
        % derivative of the Pontryaginfunction with respect to the state
        % AND costate variable(s)
        costatename=retrievediffmodelinformation(ocStruct,'costatename');
        statename=retrievediffmodelinformation(ocStruct,'statename');
        independent=retrievediffmodelinformation(ocStruct,'independent');
        statevector=cell2vectorstring([statename.value '_' independent.value costatename.value '_' independent.value 'p1']);
        pontryaginfunction=retrievediffmodelinformation(ocStruct,'specificpontryaginfunction',arcidentifier,symkernel);
        data.value=string2cell(ocmatjacobian(['[' pontryaginfunction.value ']'],statevector,symkernel),'vector');
        data.type='mathcellchar';
        data.description='derivative of the Pontryaginfunction with respect to the state and costate';

    case 'pontryaginfunctionDutp1'
        % derivative of the Pontryaginfunction with respect to the control
        % variable(s)
        if isfield(ocStruct,'pontryaginfunction') && isfield(ocStruct.pontryaginfunction,'derivative') ...
                && isfield(ocStruct.pontryaginfunction.derivative,'Dutp1')
            data.value=ocStruct.pontryaginfunction.derivative.Dutp1.term;
        else
            controlnametp1=retrievediffmodelinformation(ocStruct,'controlnametp1');
            controlvector=cell2vectorstring(controlnametp1.value);
            pontryaginfunction=retrievediffmodelinformation(ocStruct,'pontryaginfunction');
            data.value=string2cell(ocmatjacobian(['[' pontryaginfunction.value ']'],controlvector,symkernel),'vector');
        end
        data.type='mathcellchar';
        data.description='derivative of the Pontryaginfunction with respect to the control';

    case 'pontryaginfunctionDxt2'
        % second order derivative of the Pontryaginfunction with respect to
        % the state variable(s)
        if isfield(ocStruct,'pontryaginfunction') && isfield(ocStruct.pontryaginfunction,'derivative') ...
                && isfield(ocStruct.pontryaginfunction.derivative,'Dxt2')
            data.value=ocStruct.pontryaginfunction.derivative.Dxt2.term;
        else
            statenamet=retrievediffmodelinformation(ocStruct,'statenamet');
            statevector=cell2vectorstring(statenamet.value);
            pontryaginfunction=retrievediffmodelinformation(ocStruct,'pontryaginfunction');
            data.value=string2cell(ocmathessian(pontryaginfunction.value,statevector,symkernel),'matrix');
        end
        data.type='mathcellchar';
        data.description='second order derivative of the Pontryaginfunction with respect to the state';

    case 'pontryaginfunctionDutp12'
        % second order derivative of the Pontryaginfunction with respect to
        % the control variable(s)
        if isfield(ocStruct,'pontryaginfunction') && isfield(ocStruct.pontryaginfunction,'derivative') ...
                && isfield(ocStruct.pontryaginfunction.derivative,'Dutp12')
            data.value=ocStruct.pontryaginfunction.derivative.Dutp12.term;
        else
            controlnametp1=retrievediffmodelinformation(ocStruct,'controlnametp1');
            controlvector=cell2vectorstring(controlnametp1.value);
            pontryaginfunction=retrievediffmodelinformation(ocStruct,'pontryaginfunction');
            data.value=string2cell(ocmathessian(pontryaginfunction.value,controlvector,symkernel),'matrix');
        end
        data.type='mathcellchar';
        data.description='second order derivative of the Pontryaginfunction with respect to the control';

    case {'pontryaginfunctionDutp1Dxt','pontryaginfunctionDxtDutp1'}
        % second order derivative of the Pontryaginfunction with respect to
        % the control variable(s)
        if isfield(ocStruct,'pontryaginfunction') && isfield(ocStruct.pontryaginfunction,'derivative') ...
                && isfield(ocStruct.pontryaginfunction.derivative,'Dutp1Dxt')
            data.value=ocStruct.pontryaginfunction.derivative.Dutp1Dxt.term;
        else
            controlnametp1=retrievediffmodelinformation(ocStruct,'controlnametp1');
            controlvector=cell2vectorstring(controlnametp1.value);
            statenamet=retrievediffmodelinformation(ocStruct,'statenamet');
            statevector=cell2vectorstring(statenamet.value);
            pontryaginfunction=retrievediffmodelinformation(ocStruct,'pontryaginfunction');
            J=removematrixstring(ocmatjacobian(['[' pontryaginfunction.value ']'],controlvector,symkernel));
            data.value=string2cell(ocmatjacobian(J,statevector,symkernel),'matrix');
        end
        data.type='mathcellchar';
        data.description='second order derivative of the Pontryaginfunction with respect to the control';

    case 'optimalcontroldynamicsleftside'
        controlname=retrievediffmodelinformation(ocStruct,'controlname');
        controlvector=cell2vectorstring(controlname.value);
        zerolmmcindex=retrievediffmodelinformation(ocStruct,'zerolmmcindex',arcidentifier);
        lagrangemultipliercontrolname=retrievediffmodelinformation(ocStruct,'lagrangemultipliercontrolname');
        lagrangemultipliercontrolname.value(zerolmmcindex.value)=[];
        lagrangemultipliervector=cell2vectorstring(lagrangemultipliercontrolname.value);
        pontryaginfunction=retrievediffmodelinformation(ocStruct,'pontryaginfunction');
        pontryaginfunction4arc=retrievediffmodelinformation(ocStruct,'pontryaginfunction4arc',arcidentifier,symkernel);
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

    case 'invoptimalcontroldynamicsleftside'
        optimalcontroldynamicsleftside=retrievediffmodelinformation(ocStruct,'optimalcontroldynamicsleftside',arcidentifier,symkernel);
        LSM=mycell2sym(optimalcontroldynamicsleftside.value,'matrix');
        data.value=string2cell(char(inv(LSM)),'charmatrix');
        data.type='mathcellchar';
        data.description='inverse of the left side matrix of the optimal control dynamics';


    case 'optimalcontroldynamics'
        invoptimalcontroldynamicsleftside=retrievediffmodelinformation(ocStruct,'invoptimalcontroldynamicsleftside',arcidentifier,symkernel);
        pontryaginfunctionDuDX=retrievediffmodelinformation(ocStruct,'pontryaginfunctionDuDX',arcidentifier,symkernel);
        pontryaginfunctionDlmmcDX=retrievediffmodelinformation(ocStruct,'pontryaginfunctionDlmmcDX',arcidentifier,symkernel);
        implicitnonlinearcontrolindex=retrievediffmodelinformation(ocStruct,'implicitnonlinearcontrolindex',arcidentifier);
        if ~isempty(implicitnonlinearcontrolindex.value)
            canonicalsystem=retrievediffmodelinformation(ocStruct,'canonicalsystem','',symkernel);
            invLSM=mycell2sym(invoptimalcontroldynamicsleftside.value,'matrix');
            RSM1=mycell2sym(pontryaginfunctionDuDX.value,'matrix');
            if ~isempty(pontryaginfunctionDlmmcDX.value)
                RSM2=mycell2sym(pontryaginfunctionDlmmcDX.value,'matrix');
            else
                RSM2=sym([]);
            end
            RSM=[RSM1;RSM2'];
            xp1=mycell2sym(canonicalsystem.value,'vector');
            dudt=-invLSM*RSM*xp1;
            data.value=string2cell(char(dudt(implicitnonlinearcontrolindex.value)),'vector');
        else
            data.value='';
        end
        data.type='mathcellchar';
        data.description='inverse of the left side matrix of the optimal control dynamics';

    case 'tensoroptimalcontroldynamicsleftside'
        controlname=retrievediffmodelinformation(ocStruct,'controlname');
        controlvector=cell2vectorstring(controlname.value);
        zerolmmcindex=retrievediffmodelinformation(ocStruct,'zerolmmcindex',arcidentifier);
        lagrangemultipliercontrolname=retrievediffmodelinformation(ocStruct,'lagrangemultipliercontrolname');
        lagrangemultipliercontrolname.value(zerolmmcindex.value)=[];
        lagrangemultipliervector=cell2vectorstring(lagrangemultipliercontrolname.value);
        pontryaginfunction=retrievediffmodelinformation(ocStruct,'pontryaginfunction');
        pontryaginfunction4arc=retrievediffmodelinformation(ocStruct,'pontryaginfunction4arc',arcidentifier,symkernel);
        pontryaginfunctionDu2=sym(ocmathessian(pontryaginfunction.value,controlvector,symkernel));
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
        optimalvalue=retrievediffmodelinformation(ocStruct,'optimalvalue',arcidentifier);
        maximizingvariable=retrievediffmodelinformation(ocStruct,'maximizingvariable');
        equationvariablename=retrievediffmodelinformation(ocStruct,'equationvariablenameimplicit',arcidentifier);
        for ii=1:numel(maximizingvariable.value)
            optimalcontroldynamicsleftside=ocmatsubs(optimalcontroldynamicsleftside,[maximizingvariable.value{ii} '=' optimalvalue.value.(maximizingvariable.value{ii})]);
        end
        for ii=1:numel(equationvariablename.value)
            T{ii}=string2cell(removematrixstring(ocmatdiff(optimalcontroldynamicsleftside,equationvariablename.value{ii},symkernel)),'matrix');
        end
        data.value=T;
        data.type='mathcellchar';
        data.description='left side matrix of the optimal control dynamics';

    case 'jacobianoptimalcontroldynamicsrightside'
        controlname=retrievediffmodelinformation(ocStruct,'controlname');
        controlvector=cell2vectorstring(controlname.value);
        statename=retrievediffmodelinformation(ocStruct,'statename');
        costatename=retrievediffmodelinformation(ocStruct,'costatename');
        statevector=cell2vectorstring([statename.value costatename.value]);
        specificstatedynamics=retrievediffmodelinformation(ocStruct,'specificstatedynamics',arcidentifier,'maple');
        specificadjointsystem=retrievediffmodelinformation(ocStruct,'specificadjointsystem',arcidentifier,'maple');
        xp1=mycell2sym([specificstatedynamics.value specificadjointsystem.value],'vector');
        pontryaginfunction=retrievediffmodelinformation(ocStruct,'pontryaginfunction');
        J=removematrixstring(ocmatjacobian(['[' pontryaginfunction.value ']'],controlvector,symkernel));
        J=sym(ocmatjacobian(J,statevector,symkernel));
        RS=char(J*xp1);
        optimalvalue=retrievediffmodelinformation(ocStruct,'optimalvalue',arcidentifier);
        maximizingvariable=retrievediffmodelinformation(ocStruct,'maximizingvariable');
        equationvariablename=retrievediffmodelinformation(ocStruct,'equationvariablenameimplicit',arcidentifier);
        vector=cell2vectorstring(equationvariablename.value);
        for ii=1:numel(maximizingvariable.value)
            RS=ocmatsubs(RS,[maximizingvariable.value{ii} '=' optimalvalue.value.(maximizingvariable.value{ii})]);
        end
        T=string2cell(char(jacobian(sym(RS),vector)),'charmatrix');
        data.value=T;
        data.type='mathcellchar';
        data.description='left side matrix of the optimal control dynamics';

    case {'pontryaginfunctionDuDX','pontryaginfunctionDXDu'}
        try
            data.value=ocStruct.pontryaginfunction.derivative.DuDX.term;
        catch
            controlname=retrievediffmodelinformation(ocStruct,'controlname');
            controlvector=cell2vectorstring(controlname.value);
            statename=retrievediffmodelinformation(ocStruct,'statename');
            costatename=retrievediffmodelinformation(ocStruct,'costatename');
            statevector=cell2vectorstring([statename.value costatename.value]);
            pontryaginfunction=retrievediffmodelinformation(ocStruct,'pontryaginfunction');
            J=removematrixstring(ocmatjacobian(['[' pontryaginfunction.value ']'],controlvector,symkernel));
            data.value=string2cell(ocmatjacobian(J,statevector,symkernel),'matrix');
        end
        data.type='mathcellchar';
        data.description='right side of the optimal control dynamics (-Pux*dotx-Pul*dotl)';

    case {'pontryaginfunctionDuDlmmc','pontryaginfunctionDlmmcDu'}
        try
            data.value=ocStruct.pontryaginfunction.derivative.DuDlmmc.term;
        catch
            controlname=retrievediffmodelinformation(ocStruct,'controlname');
            controlvector=cell2vectorstring(controlname.value);
            lagrangemultipliercontrolname=retrievediffmodelinformation(ocStruct,'lagrangemultipliercontrolname');
            lagrangemultipliervector=cell2vectorstring(lagrangemultipliercontrolname.value);
            pontryaginfunction=retrievediffmodelinformation(ocStruct,'pontryaginfunction');
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
        statename=retrievediffmodelinformation(ocStruct,'statename');
        costatename=retrievediffmodelinformation(ocStruct,'costatename');
        statevector=cell2vectorstring([statename.value costatename.value]);
        zerolmmcindex=retrievediffmodelinformation(ocStruct,'zerolmmcindex',arcidentifier);
        lagrangemultipliercontrolname=retrievediffmodelinformation(ocStruct,'lagrangemultipliercontrolname');
        lagrangemultipliercontrolname.value(zerolmmcindex.value)=[];
        lagrangemultipliervector=cell2vectorstring(lagrangemultipliercontrolname.value);
        pontryaginfunction=retrievediffmodelinformation(ocStruct,'pontryaginfunction');
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
            controlname=retrievediffmodelinformation(ocStruct,'controlname');
            controlvector=cell2vectorstring(controlname.value);
            lagrangemultipliercontrolname=retrievediffmodelinformation(ocStruct,'lagrangemultipliercontrolname');
            lagrangemultipliervector=cell2vectorstring(lagrangemultipliercontrolname.value);
            pontryaginfunction=retrievediffmodelinformation(ocStruct,'pontryaginfunction');
            pontryaginfunctionDuDlmmc=removematrixstring(ocmatjacobian(['[' pontryaginfunction.value ']'],controlvector,symkernel));
            pontryaginfunctionDuDlmmc=ocmatjacobian(pontryaginfunctionDuDlmmc,lagrangemultipliervector,symkernel);
        end
        zerolmmcindex=retrievediffmodelinformation(ocStruct,'zerolmmcindex',arcidentifier);
        pontryaginfunctionDuDlmmc(:,zerolmmcindex.value)=[];
        data.value=string2cell(removematrixstring(char(pontryaginfunctionDuDlmmc)),'charmatrix');
        data.type='mathcellchar';
        data.description='right side of the optimal control dynamics (-Pux*dotx-Pul*dotl)';
        
    case {'zeros4optimalcontroldynamics'}
        zerolmmcindex=retrievediffmodelinformation(ocStruct,'zerolmmcindex',arcidentifier);
        lagrangemultipliercontrolname=retrievediffmodelinformation(ocStruct,'lagrangemultipliercontrolname');
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
            maximizingimplicitvariable=retrievediffmodelinformation(ocStruct,'maximizingimplicitvariable',arcidentifier);
            controlname=retrievediffmodelinformation(ocStruct,'controlname');
            zerolmmcindex=retrievediffmodelinformation(ocStruct,'zerolmmcindex',arcidentifier);
            lagrangemultipliercontrolname=retrievediffmodelinformation(ocStruct,'lagrangemultipliercontrolname');
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
        pontryaginfunction=retrievediffmodelinformation(ocStruct,'pontryaginfunction');
        zerolmmcindex=retrievediffmodelinformation(ocStruct,'zerolmmcindex',arcidentifier);
        lagrangemultipliercontrolnametp1=retrievediffmodelinformation(ocStruct,'lagrangemultipliercontrolnametp1');
        lmmcequalzero='';
        ii=0;
        for idx=zerolmmcindex.value
            ii=ii+1;
            lmmcequalzero{ii}=[lagrangemultipliercontrolnametp1.value{idx} '=0'];
        end
        zerolmmcvector=cell2vectorstring(lmmcequalzero);
        if ~isempty(zerolmmcvector)
            data.value=ocmatsubs(pontryaginfunction.value,zerolmmcvector,symkernel);
        else
            data.value=pontryaginfunction.value;
        end
        data.type='mathchar';
        data.description='specific Hamilton function for actual constraint combination';

    case 'maximizingsolution'
        % Pontryagin's maximumprinciple is used to derive explicit formulas
        % of control and Lagrange multiplier values satisfying the first
        % order necessary conditions
        solution=[];
        maximizingvariable=retrievediffmodelinformation(ocStruct,'maximizingvariable');
        maximizingexplicitvariable=retrievediffmodelinformation(ocStruct,'maximizingexplicitvariable',arcidentifier);
        maximizingimplicitvariable=retrievediffmodelinformation(ocStruct,'maximizingimplicitvariable',arcidentifier);
        zerolmmctp1=retrievediffmodelinformation(ocStruct,'zerolmmctp1',arcidentifier);
        maximizingderivativevariable=retrievediffmodelinformation(ocStruct,'maximizingderivativevariable',arcidentifier);
        controlnametp1=retrievediffmodelinformation(ocStruct,'controlnametp1');
        lagrangemultipliercontrolnametp1=retrievediffmodelinformation(ocStruct,'lagrangemultipliercontrolnametp1');
        controlnum=retrievediffmodelinformation(ocStruct,'controlnum');
        inequalitycontrolconstraintnum=retrievediffmodelinformation(ocStruct,'inequalitycontrolconstraintnum');

        pontryaginfunction4arc=retrievediffmodelinformation(ocStruct,'pontryaginfunction4arc',arcidentifier,symkernel);
        equationvariable=cell2vectorstring([maximizingexplicitvariable.value]);
        derivativevariable=cell2vectorstring([maximizingderivativevariable.value]);
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
                    solution.(maximizingimplicitvariable.value{ii})=maximizingimplicitvariable.value{ii};
                end
            end
            for jj=1:numel(solution)
                for ii=1:numel(zerolmmctp1.value)
                    % the Lagrange multipliers being zero (inactive
                    % constraint) are added
                    solution(jj).(zerolmmctp1.value{ii})='0';
                end
            end
        else
            % the implicit control values are added to the solution
            % structure, where the value is the variable itself
            for ii=1:numel(maximizingimplicitvariable.value)
                solution.(maximizingimplicitvariable.value{ii})=maximizingimplicitvariable.value{ii};
            end
            for ii=1:numel(zerolmmctp1.value)
                % the Lagrange multipliers being zero (inactive
                % constraint) are added
                solution.(zerolmmctp1.value{ii})='0';
            end
        end
        % the order of the solution fields are controlname,
        % lmmcname
        solution=orderfields(solution,maximizingvariable.value);
        for ii=1:numel(solution)
            for jj=1:controlnum.value
                solution(ii).control{jj}=solution(ii).(controlnametp1.value{jj});
            end
            for jj=1:inequalitycontrolconstraintnum.value
                solution(ii).lagrangemultcc{jj}=solution(ii).(lagrangemultipliercontrolnametp1.value{jj});
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
                && isfield(ocStruct.foc.adjointsystem.dynamics.(arc),'diff')
            data.value=ocStruct.foc.adjointsystem.dynamics.(arc).diff.term;
        else
            pontryaginfunctionDxt=retrievediffmodelinformation(ocStruct,'pontryaginfunctionDxt','',symkernel);
            statenum=retrievediffmodelinformation(ocStruct,'statenum');
            costatenamet=retrievediffmodelinformation(ocStruct,'costatenamet');
            discountratevariable=retrievediffmodelinformation(ocStruct,'discountratevariable');
            for ii=1:statenum.value
                data.value{ii}=['exp(' discountratevariable.value ')*' costatenamet.value{ii} '-(' pontryaginfunctionDxt.value{ii} ')'];
            end
        end
        data.type='mathcellchar';
        data.description='costate dynamics';

    case 'adjointsystemexplicit'
        % the adjoint system, where each cell contains one costate dynamics
        arc=arcidentifier2field(arcidentifier);
        if isfield(ocStruct,'foc') && isfield(ocStruct.foc,'adjointsystem') ...
                && isfield(ocStruct.foc.adjointsystem,'dynamics') ...
                && isfield(ocStruct.foc.adjointsystem.dynamics,arc) ...
                && isfield(ocStruct.foc.adjointsystem.dynamics.(arc),'diff')
            data.value=ocStruct.foc.adjointsystem.dynamics.(arc).diff.term;
        else
            pontryaginfunctionDxt=retrievediffmodelinformation(ocStruct,'pontryaginfunctionDxt','',symkernel);
            statenum=retrievediffmodelinformation(ocStruct,'statenum');
            costatenamet=retrievediffmodelinformation(ocStruct,'costatenamet');
            costatenametp1=retrievediffmodelinformation(ocStruct,'costatenametp1');
            discountratevariable=retrievediffmodelinformation(ocStruct,'discountratevariable');
            for ii=1:statenum.value
                dldt{ii}=['exp(' discountratevariable.value ')*' costatenamet.value{ii} '-(' pontryaginfunctionDxt.value{ii} ')'];
            end
            equationvariable=cell2vectorstring([costatenametp1.value]);
            equation=cell2vectorstring(dldt);
            % ocmatsolve is an adaptation of the native MATLAB command
            % solve for the symbolic toolbox relying
                solution=ocmatsolve(equation,equationvariable,symkernel);
            if ~isempty(solution)
                for ii=1:statenum.value
                    data.value{ii}=solution.(costatenametp1.value{ii});
                end
            else
                data.value=[];
            end
        end
        data.type='mathcellchar';
        data.description='costate dynamics';

    case 'algebraicequation'
        % if the Hamiltonian maximizing condition does not allow an
        % explicit algebraic expression for the control values, the
        % corresponding equations are added as algebaric equations yielding
        % a system of algebraic differential equations
        arc=arcidentifier2field(arcidentifier);
        if isfield(ocStruct,'foc') && isfield(ocStruct.foc,'adjointsystem') ...
                && isfield(ocStruct.foc.adjointsystem,'algebraicequation') ...
                && isfield(ocStruct.foc.adjointsystem.algebraicequation,arc)
            data.value=ocStruct.foc.adjointsystem.algebraicequation.(arc).term;
        else
            pontryaginfunction4arc=retrievediffmodelinformation(ocStruct,'pontryaginfunction4arc',arcidentifier,symkernel);
            maximizingimplicitderivativevariable=retrievediffmodelinformation(ocStruct,'maximizingimplicitderivativevariable',arcidentifier);
            if ~isempty(maximizingimplicitderivativevariable.value)
                derivativevariable=cell2vectorstring([maximizingimplicitderivativevariable.value]);
                data.value=string2cell(removematrixstring(ocmatjacobian(['[' pontryaginfunction4arc.value ']'],derivativevariable,symkernel)),'vector');
            end
        end
        data.type='mathcellchar';
        data.description='equation to determine the implicitly given control values';

    case 'specificstatedynamics'
        % the canonical system consisting of the state dynamics, adjoint
        % system and possible algebraic equations for implicitly given
        % control values.
        statedynamics=retrievediffmodelinformation(ocStruct,'statedynamics');
        optimalvalue=retrievediffmodelinformation(ocStruct,'optimalvalue',arcidentifier);
        maximizingvariable=retrievediffmodelinformation(ocStruct,'maximizingvariable');
        xp1=statedynamics.value;
        for ii=1:numel(maximizingvariable.value)
            for jj=1:numel(xp1)
                xp1{jj}=ocmatsubs(xp1{jj},[maximizingvariable.value{ii} '=' optimalvalue.value.(maximizingvariable.value{ii})]);
            end
        end
        data.value=xp1;
        data.type='mathcellchar';
        data.description='state dynamics';

    case 'specificadjointsystem'
        % the canonical system consisting of the state dynamics, adjoint
        % system and possible algebraic equations for implicitly given
        % control values.
        adjointsystem=retrievediffmodelinformation(ocStruct,'adjointsystem','',symkernel);
        optimalvalue=retrievediffmodelinformation(ocStruct,'optimalvalue',arcidentifier);
        maximizingvariable=retrievediffmodelinformation(ocStruct,'maximizingvariable');
        dldt=adjointsystem.value;

        for ii=1:numel(maximizingvariable.value)
            for jj=1:numel(dldt)
                dldt{jj}=ocmatsubs(dldt{jj},[maximizingvariable.value{ii} '=' optimalvalue.value.(maximizingvariable.value{ii})]);
            end
        end
        data.value=dldt;
        data.type='mathcellchar';
        data.description='costate dynamics';

    case 'adjointstateexplicit'
        % the canonical system consisting of the state dynamics, adjoint
        % system and possible algebraic equations for implicitly given
        % control values.
        adjointsystem=retrievediffmodelinformation(ocStruct,'adjointsystem','',symkernel);
        optimalvalue=retrievediffmodelinformation(ocStruct,'optimalvalue',arcidentifier);
        maximizingvariable=retrievediffmodelinformation(ocStruct,'maximizingvariable');
        costatenametp1=retrievediffmodelinformation(ocStruct,'costatenametp1');
        statenum=retrievediffmodelinformation(ocStruct,'statenum');
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
        maptype=retrievediffmodelinformation(ocStruct,'maptype');
        if strcmp(maptype.value,'explicit')
            data=retrievediffmodelinformation(ocStruct,'specificcanonicalsystemexplicit',arcidentifier,symkernel);
        elseif strcmp(maptype.value,'implicit')
            data=retrievediffmodelinformation(ocStruct,'specificcanonicalsystemimplicit',arcidentifier,symkernel);
        else
            data.value=[];
        end
         data.type='mathcellchar';
       
    case 'specificcanonicalsystemexplicit'
        % the canonical system consisting of the state dynamics, adjoint
        % system and possible algebraic equations for implicitly given
        % control values.
        adjointstateexplicit=retrievediffmodelinformation(ocStruct,'adjointstateexplicit',arcidentifier,symkernel);
        if isempty(adjointstateexplicit.value)
            data.value=[];
            return
        end
        statedynamicsexplicit=retrievediffmodelinformation(ocStruct,'statedynamics');
        optimalvalue=retrievediffmodelinformation(ocStruct,'optimalvalue',arcidentifier);
        maximizingvariable=retrievediffmodelinformation(ocStruct,'maximizingvariable');
        costatenametp1=retrievediffmodelinformation(ocStruct,'costatenametp1');
        statenum=retrievediffmodelinformation(ocStruct,'statenum');
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
        canonicalsystemimplicit=retrievediffmodelinformation(ocStruct,'canonicalsystemimplicit','',symkernel);
        optimalvalue=retrievediffmodelinformation(ocStruct,'optimalvalue',arcidentifier);
        maximizingvariable=retrievediffmodelinformation(ocStruct,'maximizingvariable');
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
        pontryaginfunction=retrievediffmodelinformation(ocStruct,'pontryaginfunction','',symkernel);
        optimalvalue=retrievediffmodelinformation(ocStruct,'optimalvalue',arcidentifier,symkernel);
        maximizingvariable=retrievediffmodelinformation(ocStruct,'maximizingvariable');
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
        adjointsystem=retrievediffmodelinformation(ocStruct,'adjointsystemexplicit',arcidentifier,symkernel);
        statedynamics=retrievediffmodelinformation(ocStruct,'statedynamics');
        xp1=[statedynamics.value adjointsystem.value];
        data.value=xp1;
        data.type='mathcellchar';
        data.description='state dynamics, adjoint system';

    case 'canonicalsystemimplicit'
        % the canonical system consisting of the state dynamics, adjoint
        % system and possible algebraic equations for implicitly given
        % control values.
        adjointsystem=retrievediffmodelinformation(ocStruct,'adjointsystem','',symkernel);
        statedynamicsimplicit=retrievediffmodelinformation(ocStruct,'statedynamicsimplicit');
        xp1=[statedynamicsimplicit.value adjointsystem.value];
        data.value=xp1;
        data.type='mathcellchar';
        data.description='state dynamics, adjoint system';

    case 'canonicalsystem'
        % the canonical system consisting of the state dynamics, adjoint
        % system and possible algebraic equations for implicitly given
        % control values.
        adjointsystem=retrievediffmodelinformation(ocStruct,'adjointsystem','',symkernel);
        statedynamics=retrievediffmodelinformation(ocStruct,'statedynamics');
        xp1=[statedynamics.value adjointsystem.value];
        data.value=xp1;
        data.type='mathcellchar';
        data.description='state dynamics, adjoint system';

    case 'generalcanonicalsystem'
        % the canonical system consisting of the state dynamics, adjoint
        % system and possible algebraic equations for implicitly given
        % control values.
        adjointsystem=retrievediffmodelinformation(ocStruct,'adjointsystem','',symkernel);
        algebraicequation=retrievediffmodelinformation(ocStruct,'algebraicequation',arcidentifier,symkernel);
        statedynamics=retrievediffmodelinformation(ocStruct,'statedynamics');
        xp1=[statedynamics.value adjointsystem.value algebraicequation.value];
        data.value=xp1;
        data.type='mathcellchar';
        data.description='state dynamics, adjoint system (and algebraic equations for implicit control values)';

    case 'canonicalsystemjacobian'
        % Jacobian of the canonical system with respect to state, costate
        % and implicit control/Lagrange multiplier variables
        specificcanonicalsystem=retrievediffmodelinformation(ocStruct,'specificcanonicalsystemimplicit',arcidentifier,symkernel);
        data.description='Jacobian of the implicit canonical system map with respect to state, costate and implicit control/Lagrange multiplier variables';
        canonicalsystem=cell2vectorstring(specificcanonicalsystem.value);
        equationvariablenamet=retrievediffmodelinformation(ocStruct,'equationvariablenamet',arcidentifier);
        equationvariablenametp1=retrievediffmodelinformation(ocStruct,'equationvariablenametp1',arcidentifier);
        variable=cell2vectorstring([equationvariablenamet.value equationvariablenametp1.value]);
        data.value=string2cell(removematrixstring(ocmatjacobian(canonicalsystem,variable,symkernel)),'matrix');
        data.type='mathcellchar';


    case 'statecostatejacobian'
        % Jacobian of the canonical system with respect to state, costate
        % and implicit control/Lagrange multiplier variables
        specificstatedynamics=retrievediffmodelinformation(ocStruct,'specificstatedynamics',arcidentifier,symkernel);
        specificadjointsystem=retrievediffmodelinformation(ocStruct,'specificadjointsystem',arcidentifier,symkernel);
        xp1=cell2vectorstring([specificstatedynamics.value specificadjointsystem.value]);
        equationvariablename=retrievediffmodelinformation(ocStruct,'equationvariablename',arcidentifier);
        variable=cell2vectorstring(equationvariablename.value);
        data.value=string2cell(removematrixstring(ocmatjacobian(xp1,variable,symkernel)),'matrix');
        data.type='mathcellchar';
        data.description='Jacobian of the canonical system with respect to state, costate and implicit control/Lagrange multiplier variables';

    case 'specificstatedynamicsjacobian'
        % Jacobian of the canonical system with respect to state, costate
        % and implicit control/Lagrange multiplier variables
        specificstatedynamics=retrievediffmodelinformation(ocStruct,'specificstatedynamics',arcidentifier,symkernel);
        xp1=cell2vectorstring(specificstatedynamics.value);
        equationvariablenamet=retrievediffmodelinformation(ocStruct,'equationvariablenamet',arcidentifier);
        equationvariablenametp1=retrievediffmodelinformation(ocStruct,'equationvariablenametp1',arcidentifier);
        variable=cell2vectorstring([equationvariablenamet.value equationvariablenametp1.value]);
        data.value=string2cell(removematrixstring(ocmatjacobian(xp1,variable,symkernel)),'matrix');
        data.type='mathcellchar';
        data.description='Jacobian of the canonical system with respect to state, costate and implicit control/Lagrange multiplier variables';

    case 'canonicalsystemhessian'
        % Jacobian of the canonical system with respect to state, costate
        % and implicit control/Lagrange multiplier variables
        canonicalsystem=retrievediffmodelinformation(ocStruct,'specificcanonicalsystem',arcidentifier,symkernel);
        equationvariablename=retrievediffmodelinformation(ocStruct,'equationvariablename',arcidentifier);
        variable=cell2vectorstring(equationvariablename.value);
        xp1=canonicalsystem.value;
        for ii=1:numel(xp1)
            H{ii}=string2cell(removematrixstring(ocmathessian(xp1{ii},variable,symkernel)),'matrix');
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
        maptype=retrievediffmodelinformation(ocStruct,'maptype');
        if strcmp(maptype.value,'explicit')
            specificcanonicalsystem=retrievediffmodelinformation(ocStruct,'specificcanonicalsystemexplicit',arcidentifier,symkernel);
            data.description='Jacobian of the explicit canonical system map with respect to exogenous parameter variables';
        elseif strcmp(maptype.value,'implicit')
            specificcanonicalsystem=retrievediffmodelinformation(ocStruct,'specificcanonicalsystemimplicit',arcidentifier,symkernel);
            data.description='Jacobian of the implicit canonical system map with respect to exogenous parameter variables';
        else
            data.value=[];
        end
        xp1=cell2vectorstring(specificcanonicalsystem.value);
        parametername=retrievediffmodelinformation(ocStruct,'parametername',arcidentifier);
        variable=cell2vectorstring(parametername.value);
        data.value=string2cell(removematrixstring(ocmatjacobian(xp1,variable,symkernel)),'matrix');
        data.type='mathcellchar';
        data.description='Jacobian of the canonical system with respect to exogenous parameter variables';

    case 'canonicalsystemderivativetime'
        % Jacobian of the canonical system with respect to exogenous
        % parameter variables
        canonicalsystem=retrievediffmodelinformation(ocStruct,'specificcanonicalsystem',arcidentifier,symkernel);
        xp1=cell2vectorstring(canonicalsystem.value);
        timevariable=retrievediffmodelinformation(ocStruct,'independent');
        variable=cell2vectorstring(timevariable.value);
        data.value=string2cell(removematrixstring(ocmatjacobian(xp1,variable,symkernel)),'matrix');
        data.type='mathcellchar';
        data.description='Derivative of the canonical system with respect to the independent variable';

    case 'canonicalsystemtotalhessian'
        % Jacobian of the canonical system with respect to exogenous
        % parameter variables
        canonicalsystem=retrievediffmodelinformation(ocStruct,'specificcanonicalsystem',arcidentifier,symkernel);
        equationvariablename=retrievediffmodelinformation(ocStruct,'equationvariablename',arcidentifier);
        parametername=retrievediffmodelinformation(ocStruct,'parametername',arcidentifier);
        variable=cell2vectorstring([equationvariablename.value(:);parametername.value]);
        xp1=canonicalsystem.value;
        for ii=1:numel(xp1)
            H{ii}=string2cell(removematrixstring(ocmathessian(xp1{ii},variable,symkernel)),'matrix');
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
        controlnametp1=retrievediffmodelinformation(ocStruct,'controlnametp1');
        controlnum=retrievediffmodelinformation(ocStruct,'controlnum');
        lagrangemultipliercontrolnametp1=retrievediffmodelinformation(ocStruct,'lagrangemultipliercontrolnametp1');
        inequalitycontrolconstraintnum=retrievediffmodelinformation(ocStruct,'inequalitycontrolconstraintnum');

        for ii=1:controlnum.value
            value.(controlnametp1.value{ii})=ocStruct.foc.value.control.(arcidentifier2field(arcidentifier)).term{ii};
        end
        for ii=1:inequalitycontrolconstraintnum.value
            value.(lagrangemultipliercontrolnametp1.value{ii})=ocStruct.foc.value.lagrangemultcc.(arcidentifier2field(arcidentifier)).term{ii};
        end
        data.type='structmathchar';
        data.value=value;
        data.description='optimal control values as derived from the Hamiltonian maximizing condition';

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

