function b=standardmodeltestconsistency(ocStruct,testtype,symkernel,sectionStruct,sectiontext)
%
% returns 0 if input is consistent, 1 otherwise
b=0;

switch testtype
    case 'maininitfile'
        [mainpropnames,mandatorysectionname]=standardmodelproperties('init');
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
        % test if 'control' appears in variable section
        o=regexpi(sectiontext(sectionStruct.variable(1):sectionStruct.variable(2)),'\<control\:{2}');
        if isempty([o{:}])
            b=1;
            ocmatmsg('''control'' is missing in section variable.');
            return
        end
        % test if 'ode' appears in statedynamics section
        o=regexpi(sectiontext(sectionStruct.statedynamics(1):sectionStruct.statedynamics(2)),'\<ode\:{2}');
        if isempty([o{:}])
            b=1;
            ocmatmsg('''ode'' is missing in section Statedynamics.');
            return
        end
        % test if 'int' appears in objective section
        o=regexpi(sectiontext(sectionStruct.objective(1):sectionStruct.objective(2)),'\<int\:{2}');
        if isempty([o{:}])
            b=1;
            ocmatmsg('''int'' is missing in section Objective.');
            return
        end
        
    case 'focinitfile'
        [mainpropnames mandatorysectionname]=standardmodelproperties('foc');
        % test if mandatory sections exist and if these sections are not
        % empty
        for name=mandatorysectionname
            if ~isfield(sectionStruct,name{1}) || diff(sectionStruct.(name{1}))<0
                b=1;
                ocmatmsg(['Section ''' name{1} ''' is missing. Or syntactically incorrect.']);
                return
            end
        end
        
    case 'basicmodelstructure'
        % consistency test after succesfully processing the initialization
        % file
        
        % test if number of states and costates are equal
        arcnum=retrievemodelinformation(ocStruct,'arcnum');
        arcidentifier=retrievemodelinformation(ocStruct,'arcidentifier');
        statenum=retrievemodelinformation(ocStruct,'statenum');
        statename=retrievemodelinformation(ocStruct,'statename');
        costatenum=retrievemodelinformation(ocStruct,'costatenum');
        costatename=retrievemodelinformation(ocStruct,'costatename');
        variationstatename=retrievemodelinformation(ocStruct,'variationstatename');
        variationstatenum=retrievemodelinformation(ocStruct,'variationstatenum');
        variationcostatename=retrievemodelinformation(ocStruct,'variationcostatename');
        variationcostatenum=retrievemodelinformation(ocStruct,'variationcostatenum');
        variationcontrolname=retrievemodelinformation(ocStruct,'variationcontrolname');
        variationcontrolnum=retrievemodelinformation(ocStruct,'variationcontrolnum');
        variationparameternum=retrievemodelinformation(ocStruct,'variationparameternum');
        controlnum=retrievemodelinformation(ocStruct,'controlnum');
        controlname=retrievemodelinformation(ocStruct,'controlname');
        parameternum=retrievemodelinformation(ocStruct,'parameternum');
        parametername=retrievemodelinformation(ocStruct,'parametername');
        parametervalue=retrievemodelinformation(ocStruct,'parametervalue');
        independent=retrievemodelinformation(ocStruct,'independent');
        endtime=retrievemodelinformation(ocStruct,'endtime');
        statedynamics=retrievemodelinformation(ocStruct,'statedynamics');
        discountratevariable=retrievemodelinformation(ocStruct,'discountratevariable');
        variablename=retrievemodelinformation(ocStruct,'variablename');
        stdvariablename=standardmodelvariables();
        
        objectiveintegrand=retrievemodelinformation(ocStruct,'objectiveintegrand');
        salvagevalue=retrievemodelinformation(ocStruct,'salvagevalue');
        endtimediscountfactor=retrievemodelinformation(ocStruct,'endtimediscountfactor');
        discountfactor=retrievemodelinformation(ocStruct,'discountfactor');
        endtimediscountratevariable=retrievemodelinformation(ocStruct,'endtimediscountratevariable');
        inequalitycontrolconstraintidentifier=retrievemodelinformation(ocStruct,'inequalitycontrolconstraintidentifier');
        inequalitycontrolconstraint=retrievemodelinformation(ocStruct,'inequalitycontrolconstraint');
        inequalitycontrolconstraintnum=retrievemodelinformation(ocStruct,'inequalitycontrolconstraintnum');
        lagrangemultipliercontrolname=retrievemodelinformation(ocStruct,'lagrangemultipliercontrolname');
        inequalitystateconstraintidentifier=retrievemodelinformation(ocStruct,'inequalitystateconstraintidentifier');
        inequalitystateconstraint=retrievemodelinformation(ocStruct,'inequalitystateconstraint');
        inequalitystateconstraintnum=retrievemodelinformation(ocStruct,'inequalitystateconstraintnum');
        lagrangemultiplierstatename=retrievemodelinformation(ocStruct,'lagrangemultiplierstatename');
        exogenousfunctionnum=retrievemodelinformation(ocStruct,'exogenousfunctionnum');
        exogenousfunctionname=retrievemodelinformation(ocStruct,'exogenousfunctionname');
        exogenousfunctionterm=retrievemodelinformation(ocStruct,'exogenousfunctionterm');
        objectivetype=retrievemodelinformation(ocStruct,'objectivetype');
        parametervalueexpression=retrievemodelinformation(ocStruct,'parametervalueexpression');
        dependentvar=getbasicname('dependent');
        parametervar=getbasicname('parametervariables');
        arcidentifiervar=getbasicname('arcidentifiervar');
        
        % test if the symbolic expressions are syntactically correct
        testsym=testsymbolic(statedynamics.value,symkernel);
        if ~all(testsym)
            b=1;
            ocmatmsg('Line %d of the state dynamics is syntactically incorrect.\n',find(~testsym));
            return
        end
        if any(strcmp(objectivetype.value,'integral'))
            testsym=testsymbolic(objectiveintegrand.value,symkernel);
            if ~(testsym)
                b=1;
                ocmatmsg('Objective integrand is syntactically incorrect.\n');
                return
            end
        end
        if ~isempty(discountfactor.value)
            if ~testsymbolic(discountfactor.value,symkernel)
                b=1;
                ocmatmsg('Discount factor for integrand is syntactically incorrect.\n');
                return
            end
        end
        if inequalitycontrolconstraintnum.value
            testsym=testsymbolic(inequalitycontrolconstraint.value,symkernel);
            if ~all(testsym)
                b=1;
                ocmatmsg('Inequality control contraint %d is syntactically incorrect.\n',find(~testsym));
                return
            end
        end
        if inequalitystateconstraintnum.value
            testsym=testsymbolic(inequalitystateconstraint.value,symkernel);
            if ~all(testsym)
                b=1;
                ocmatmsg('Inequality state contraint %d is syntactically incorrect.\n',find(~testsym));
                return
            end
        end
        if ~isempty(salvagevalue.value)
            if ~testsymbolic(salvagevalue.value,symkernel)
                b=1;
                ocmatmsg('Salvagevalue function is syntactically incorrect.\n');
                return
            end
        end
        if ~isempty(endtimediscountfactor.value)
            if ~testsymbolic(endtimediscountfactor.value,symkernel)
                b=1;
                ocmatmsg('Discount factor for salvagevalue is syntactically incorrect.\n');
                return
            end
        end
        if statenum.value~=costatenum.value
            b=1;
            ocmatmsg('Number of states %d and costates %d are different.\n',statenum.value,costatenum.value);
            return
        end
        % test if number of states and number of ODEs are equal
        if statenum.value~=numel(statedynamics.value)
            b=1;
            ocmatmsg('Number of states %d and number of ODEs of the statedynamics %d are different.\n',statenum.value,numel(statedynamics.value));
            return
        end
        
        if ~isempty(salvagevalue.value)
            symbols=ocmatfindsym([statedynamics.value objectiveintegrand.value inequalitycontrolconstraint.value inequalitystateconstraint.value {salvagevalue.value} exogenousfunctionterm.value parametervalueexpression.value],symkernel);
        else
            symbols=ocmatfindsym([statedynamics.value objectiveintegrand.value inequalitycontrolconstraint.value inequalitystateconstraint.value exogenousfunctionterm.value parametervalueexpression.value],symkernel);
        end
        
        for ii=1:exogenousfunctionnum.value
            symbols(strcmp(symbols,exogenousfunctionname.value{ii}))=[];
        end
        if variationparameternum.value
            for ii=1:variationstatenum.value
                symbols(strcmp(symbols,variationstatename.value{ii}))=[];
            end
            for ii=1:variationcostatenum.value
                symbols(strcmp(symbols,variationcostatename.value{ii}))=[];
            end
            for ii=1:variationcontrolnum.value
                symbols(strcmp(symbols,variationcontrolname.value{ii}))=[];
            end
        end
        for ii=1:controlnum.value
            fidx=strcmp(symbols,controlname.value{ii});
            if ~any(fidx)
                b=1;
                ocmatmsg('Control ''%s'' does not appear in the model functions.\n',controlname.value{ii});
                return
            end
            symbols(fidx)=[];
        end
%         for ii=1:statenum.value
%             fidx=strcmp(symbols,statename.value{ii});
%             if ~any(fidx)
%                 ocmatmsg('State ''%s'' does not appear in the model functions.\n',statename.value{ii});
%                 if interruptquestion()
%                     b=2;
%                     return
%                 end
%             end
%             symbols(fidx)=[];
%         end

        fn=fieldnames(ocStruct.variable);
        for ii=1:length(fn)
            varname=ocStruct.variable.(fn{ii});
            if iscell(varname.name)
                for jj=1:length(varname.name)
                    fidx=strcmp(symbols,varname.name{jj});
                    symbols(fidx)=[];
                end
            else
                fidx=strcmp(symbols,varname.name);
                symbols(fidx)=[];
            end
        end
        % remove independent variable from list of symbols
        symbols(strcmp(symbols,independent.value))=[];
        % remove endtime variable from list of symbols
        symbols(strcmp(symbols,endtime.value))=[];

        % check if name of costate appear as a symbol
        for ii=1:costatenum.value
            if any(strcmp(symbols,costatename.value{ii}))
                b=1;
                ocmatmsg('Costate name ''%s'' appears as a symbol in the model functions.\n',costatename.value{ii});
                return
            end
        end
        
        % check if name of the dependent variable appear as a symbol
        internalvar={dependentvar,arcidentifiervar,parametervar};
        internalvartype={'dependent variable','variable for arcidentifer','variable for parameter vecttor'};
        for ii=1:numel(internalvar)
            if any(strcmp(symbols,internalvar{ii}))
                b=1;
                ocmatmsg('%s ''%s'' appears as a symbol in the model functions.\n',internalvartype{ii},internalvar{ii});
                return
            end
        end

        % warning if discountrate variable is used outside the discount
        % factor
        %         if any(strcmp(symbols,discountratevariable.value))
        %             ocmatmsg('The variable for the discount rate ''%s'' appears as a symbol in the model functions.\n',discountratevariable.value);
        %             if interruptquestion()
        %                 b=2;
        %                 return
        %             end
        %         end

        % check if Lagrange multiplier variable appears as a symbol
        for ii=1:inequalitycontrolconstraintnum.value
            if any(strcmp(symbols,lagrangemultipliercontrolname.value{ii}))
                b=1;
                ocmatmsg('The Lagrange multiplier variable ''%s'' appears as a symbol in the model functions.\n',internalvartype{ii},internalvar{ii});
                return
            end
        end

        % check if Lagrange multiplier variable appears as a symbol
        for ii=1:inequalitystateconstraintnum.value
            if any(strcmp(symbols,lagrangemultiplierstatename.value{ii}))
                b=1;
                ocmatmsg('The Lagrange multiplier variable ''%s'' appears as a symbol in the model functions.\n',internalvartype{ii},internalvar{ii});
                return
            end
        end
        
        % check if every symbol appears in the parameter list of the
        % initialization file
        symbols=[symbols ocmatfindsym(discountratevariable.value,symkernel)];
        if ~isempty(endtimediscountratevariable.value) && ~strcmp(discountratevariable.value,endtimediscountratevariable.value)
            symbols=[symbols ocmatfindsym(endtimediscountratevariable.value,symkernel)];
        end
        for ii=1:1:numel(symbols)
            if ~any(strcmp(parametername.value,symbols{ii}))
                b=1;
                ocmatmsg('Parameter ''%s'' does not appear in the parameter list.\n',symbols{ii});
                return
            end
        end

        % check that objective value function is not empty
        if any(strcmp(objectivetype.value,'integral'))
            if isempty(objectiveintegrand.value)
                b=1;
                ocmatmsg('Integrand of the objective value function is empty.\n');
                return
            end
        end
        
        % return message if number of symbols and number of parameter
        % variables are not equal
%         if numel(symbols)~=parameternum.value
%             for ii=1:numel(symbols)
%                 parametername.value(strcmp(parametername.value,symbols{ii}))=[];
%             end
%             ocmatmsg('Number of detected symbols is smaller than number of parameter in the list.\nParameter:\n')
%             ocmatmsg('''%s''\n',parametername.value{:});
%             ocmatmsg('in initialization file may never be used.\n')
%             if interruptquestion()
%                 b=2;
%                 return
%             end
%         end

        % test if provided parameter values are finite
        fidx=isnan(parametervalue.value);
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
            ocmatmsg('For a list of the standard variable names type ''standardmodelvariables()''\n');
            if interruptquestion()
                b=2;
                return
            end
        end

        % test that constraint identifiers appearing the constraint
        % combinations refer to an actual constraint identifier
        if inequalitycontrolconstraintnum.value
            for ii=1:arcnum.value
                ccomb=retrievemodelinformation(ocStruct,'constraintcombination',arcidentifier.value{ii});
                % make cells
                ccomb.value=regexp(ccomb.value,'[\ ]*','split');
                for jj=1:numel(ccomb.value)
                    o=regexp(inequalitycontrolconstraintidentifier.value,['\<' ccomb.value{jj}{1} '\>']);
                    if isempty(regexp(ccomb.value{jj},'\<\[\]','ONCE')) && ...
                            isempty([o{:}])
                        b=1;
                        ocmatmsg('Identifier ''%s'' for the constraint combiantion of arc ''%s'' is not a defined constraint identifier.\n',ccomb.value{jj},arcidentifier.value{ii});
                        return
                    end
                end
            end
        end

        % test that constraint identifiers appearing the constraint
        % combinations refer to an actual constraint identifier
        if inequalitystateconstraintnum.value
            for ii=1:arcnum.value
                scomb=retrievemodelinformation(ocStruct,'constraintcombination',arcidentifier.value{ii});
                % make cells
                scomb.value=regexp(scomb.value,'[\ ]*','split');
                for jj=1:numel(scomb.value)
                    o=regexp(inequalitystateconstraintidentifier.value,['\<' scomb.value{jj}{1} '\>']);
                    if isempty(regexp(scomb.value{jj},'\<\[\]','ONCE')) && ...
                            isempty([o{:}])
                        b=1;
                        ocmatmsg('Identifier ''%s'' for the constraint combiantion of arc ''%s'' is not a defined constraint identifier.\n',scomb.value{jj},arcidentifier.value{ii});
                        return
                    end
                end
            end
        end
end

