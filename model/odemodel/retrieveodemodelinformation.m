function data=retrieveodemodelinformation(ocStruct,propertyclass,arcidentifier,symkernel)
%
% RETRIEVEMODELINFORMATION returns a structure with information about the
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

    case 'parameternum'
        % number of parameter(s)
        data.type='integer';
        if isfield(ocStruct,'parameter')
            data.value=ocStruct.parameter.num;
        else
            data.value=0;
        end
        data.description='number of exogenous parameters';

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
            exogenousfunctionname=retrieveodemodelinformation(ocStruct,'exogenousfunctionname');
            exogenousfunctionargument=retrieveodemodelinformation(ocStruct,'exogenousfunctionargument');
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
        exogenousfunctionname=retrieveodemodelinformation(ocStruct,'exogenousfunctionname');
        if ~isempty(exogenousfunctionname.value)
            for ii=1:length(exogenousfunctionname.value)
                data.value{ii}=ocStruct.exogenousfunction.(exogenousfunctionname.value{ii}).argument;
            end
        else
            data.value=[];
        end
        data.description='arguments of exogenous functions';
    case 'exogenousfunctionterm'
        % variable name of the control(s) as defined by the user in the
        % initialization file
        data.type='cellchar';
        exogenousfunctionname=retrieveodemodelinformation(ocStruct,'exogenousfunctionname');
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
        parametername=retrieveodemodelinformation(ocStruct,'parametername');
        data.value=zeros(1,numel(parametername.value));
        for ii=1:numel(parametername.value)
            data.value(ii)=ocStruct.parameter.variable.(parametername.value{ii});
        end
        data.description='vector of user provided parameter values';

    case 'algebraicequationnum'
        data.type='integer';
        data.value=0;
        data.description='number of algebraic equations';

    case 'variablename'
        % general variable name(s) of the optimal control problem
        data.type='cellchar';
        data.value=fieldnames(ocStruct.variable);
        data.description='general variable name';

    case 'equationvariablename'
        % dependent variables of the canonical system for a specific arc
        data.type='cellchar';
        statename=retrievemodelinformation(ocStruct,'statename');
        data.value=statename.value;
        data.description='dependent variables of the canonical system';

    case 'statenum'
        % the number of state(s)
        data.type='integer';
        data.value=ocStruct.variable.state.num;
        data.description='number of state';

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
        data.value=ocStruct.variable.state.num;
        data.description='number of ODEs for the canonical system';

    case 'dynamics'
        % ODEs describing the evolution of the state(s)
        data.type='mathchar';
        data.value=ocStruct.constraint.derivative.state.term;
        data.description='ODEs of the dynamics';

    case 'dynamicsjacobian'
        % Jacobian of the canonical system with respect to state, costate
        % and implicit control/Lagrange multiplier variables
        dynamics=retrieveodemodelinformation(ocStruct,'dynamics');
        dxdt=cell2vectorstring(dynamics.value);
        statename=retrieveodemodelinformation(ocStruct,'statename');
        variable=cell2vectorstring(statename.value);
        data.value=string2cell(removematrixstring(ocmatjacobian(dxdt,variable,symkernel)),'matrix');
        data.type='mathcellchar';
        data.description='Jacobian of the dynamics with respect to state';

    case 'dynamicsparameterjacobian'
        % Jacobian of the canonical system with respect to state, costate
        % and implicit control/Lagrange multiplier variables
        dynamics=retrieveodemodelinformation(ocStruct,'dynamics');
        dxdt=cell2vectorstring(dynamics.value);
        parametername=retrieveodemodelinformation(ocStruct,'parametername');
        variable=cell2vectorstring(parametername.value);
        data.value=string2cell(removematrixstring(ocmatjacobian(dxdt,variable,symkernel)),'matrix');
        data.type='mathcellchar';
        data.description='Jacobian of the dynamics with respect to parametervalues';
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

