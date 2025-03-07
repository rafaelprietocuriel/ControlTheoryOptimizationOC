function b=odemodeltestconsistency(ocStruct,testtype,symkernel,sectionStruct,sectiontext)
%
% returns o if input is consistent, 1 otherwise
b=0;

switch testtype
    case 'maininitfile'
        [mainpropnames mandatorysectionname]=odemodelproperties('init');
        % consistency test for the syntactically correctness of the
        % initialization file
        %
        
        % test if mandatory sections exist and if these sections are not
        % empty
        for name=mandatorysectionname
            if ~isfield(sectionStruct,name{1}) || diff(sectionStruct.(name{1}))<0
                b=1;
                ocmatmsg(['Section ''' name{1} ''' is missing. Or syntactically incorrect.']);
                return
            end
        end
        % test if 'state' appears in variable section
        o=regexpi(sectiontext(sectionStruct.variable(1):sectionStruct.variable(2)),'\<state\:{2}');
        if isempty([o{:}])
            b=1;
            ocmatmsg('''state'' is missing in section variable.');
            return
        end
        % test if 'ode' appears in dynamics section
        o=regexpi(sectiontext(sectionStruct.dynamics(1):sectionStruct.dynamics(2)),'\<ode\:{2}');
        if isempty([o{:}])
            b=1;
            ocmatmsg('''ode'' is missing in section Statedynamics.');
            return
        end
        
        
    case 'basicmodelstructure'
        % consistency test after succesfully processing the initialization
        % file
        
        % test if number of states and costates are equal
        statenum=retrieveodemodelinformation(ocStruct,'statenum');
        statename=retrieveodemodelinformation(ocStruct,'statename');
        parameternum=retrieveodemodelinformation(ocStruct,'parameternum');
        parametername=retrieveodemodelinformation(ocStruct,'parametername');
        parametervalue=retrieveodemodelinformation(ocStruct,'parametervalue');
        independent=retrieveodemodelinformation(ocStruct,'independent');
        dynamics=retrieveodemodelinformation(ocStruct,'dynamics');
        variablename=retrieveodemodelinformation(ocStruct,'variablename');
        stdvariablename=odemodelvariables();
        
        exogenousfunctionnum=retrieveodemodelinformation(ocStruct,'exogenousfunctionnum');
        exogenousfunctionname=retrieveodemodelinformation(ocStruct,'exogenousfunctionname');
        exogenousfunctionterm=retrieveodemodelinformation(ocStruct,'exogenousfunctionterm');
        dependentvar=getbasicname('dependent');
        parametervar=getbasicname('parametervariables');
        
        % test if the symbolic expressions are syntactically correct
        testsym=testsymbolic(dynamics.value,symkernel);
        if ~all(testsym)
            b=1;
            ocmatmsg('Line %d of the state dynamics is syntactically incorrect.\n',find(~testsym));
            return
        end
        % test if number of states and number of ODEs are equal
        if statenum.value~=numel(dynamics.value)
            b=1;
            ocmatmsg('Number of states %d and number of ODEs of the dynamics %d are different.\n',statenum.value,numel(dynamics.value));
            return
        end
        symbols=ocmatfindsym([dynamics.value exogenousfunctionterm.value],symkernel);

        for ii=1:exogenousfunctionnum.value
            symbols(strcmp(symbols,exogenousfunctionname.value{ii}))=[];
        end
        for ii=1:statenum.value
            fidx=strcmp(symbols,statename.value{ii});
            if ~any(fidx)
                ocmatmsg('State ''%s'' does not appear in the model functions.\n',statename.value{ii});
                if interruptquestion()
                    b=2;
                    return
                end
            end
            symbols(fidx)=[];
        end
        % remove independent variable from list of symbols
        symbols(strcmp(symbols,independent.value))=[];
        
        % check if name of the dependent variable appear as a symbol
        internalvar={dependentvar,parametervar};
        internalvartype={'dependent variable','variable for parameter vecttor'};
        for ii=1:numel(internalvar)
            if any(strcmp(symbols,internalvar{ii}))
                b=1;
                ocmatmsg('%s ''%s'' appears as a symbol in the model functions.\n',internalvartype{ii},internalvar{ii});
                return
            end
        end
        
        for ii=1:1:numel(symbols)
            if ~any(strcmp(parametername.value,symbols{ii}))
                b=1;
                ocmatmsg('Parameter ''%s'' does not appear in the parameter list.\n',symbols{ii});
                return
            end
        end
        
        % return message if number of symbols and number of parameter
        % variables are not equal
        if numel(symbols)~=parameternum.value
            for ii=1:numel(symbols)
                parametername.value(strcmp(parametername.value,symbols{ii}))=[];
            end
            ocmatmsg('Number of detected symbols is smaller than number of parameter in the list.\nParameter:\n')
            ocmatmsg('''%s''\n',parametername.value{:});
            ocmatmsg('in initialization file may never be used.\n')
            if interruptquestion()
                b=2;
                return
            end
        end

        % test if provided parameter values are finite
        fidx=~isfinite(parametervalue.value);
        if any(fidx)
            fidx=find(fidx);
            for ii=fidx
                ocmatmsg('Parameter value of %s is not a finite number\n',parametername.value{ii});
                ocmatmsg('Possibly no value provided in the initialization file\n');
            end
            if interruptquestion()
                b=2;
                return
            end
        end

        % test if variable names defined in the initialization file are
        % standard variable names for the optimal control problem
        for ii=numel(variablename.value):-1:1
            if any(strcmp(stdvariablename,variablename.value{ii}))
                variablename.value(ii)=[];
            end
        end
        if ~isempty(variablename.value)
            ocmatmsg('Variable name ''%s'' is not a standard variable name. Maybe misspelled?\n',variablename.value{:})
            ocmatmsg('For a list of the standard variable names type ''odemodelvariables()''\n');
            if interruptquestion()
                b=2;
                return
            end
        end
end

