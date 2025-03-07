function ocStruct=standardmodelsupplement(ocStruct,supplementtype)
%
%
switch supplementtype
    case 'maininitfile'
        % if no arc is defined set unconstrained case
        if ~isfield(ocStruct,'arc') || ~isfield(ocStruct.arc,'identifier') ...
                || isempty(ocStruct.arc.identifier)
            ocStruct.arc.identifier{1}='0';
            ocStruct.arc.argument=0;
            ocStruct.arc.constraintcombination{1}='[]';
            ocStruct.arc.controlvaluecombination{1}='1';
            ocStruct.arc.num=1;
        end

        % set default variable name for independent variable
        if ~isfield(ocStruct.variable,'independent') || ~isfield(ocStruct.variable.independent,'name') ...
                || isempty(ocStruct.variable.independent.name)
            ocStruct.variable.independent.name=getbasicname('independent');
            ocStruct.variable.independent.num=1;
        end
        % set default costate name if not set by the user
        if ~isfield(ocStruct.variable,'costate') || ~isfield(ocStruct.variable.costate,'name') ...
                || isempty(ocStruct.variable.costate.name)
            statenum=retrievemodelinformation(ocStruct,'statenum');
            basename=getbasicname('costate');
            for ii=1:statenum.value
                ocStruct.variable.costate.name{ii}=[basename num2str(ii)];
            end
            ocStruct.variable.costate.num=statenum.value;
        end
        % set default variable name for end time variable
        if ~isfield(ocStruct.variable,'endtime') || ~isfield(ocStruct.variable.endtime,'name') ...
                || isempty(ocStruct.variable.endtime.name)
            ocStruct.variable.endtime.name=getbasicname('endtime');
            ocStruct.variable.endtime.num=1;
        end

        % set default variable name for Lagrange multiplier if control constraints
        % are given
        inequalitycontrolconstraintnum=retrievemodelinformation(ocStruct,'inequalitycontrolconstraintnum');
        if inequalitycontrolconstraintnum.value
            if ~isfield(ocStruct.variable,'lagrangemultcc') || ~isfield(ocStruct.variable.lagrangemultcc,'name') ...
                    || isempty(ocStruct.variable.lagrangemultcc.name)
                basename=getbasicname('lagrangemultcc');
                for ii=1:inequalitycontrolconstraintnum.value
                    ocStruct.variable.lagrangemultcc.name{ii}=[basename num2str(ii)];
                end
                ocStruct.variable.lagrangemultcc.num=inequalitycontrolconstraintnum.value;
            end
        end

        % set default variable name for Lagrange multiplier if state constraints
        % are given
        inequalitystateconstraintnum=retrievemodelinformation(ocStruct,'inequalitystateconstraintnum');
        if inequalitystateconstraintnum.value
            if ~isfield(ocStruct.variable,'lagrangemultsc') || ~isfield(ocStruct.variable.lagrangemultsc,'name') ...
                    || isempty(ocStruct.variable.lagrangemultsc.name)
                basename=getbasicname('lagrangemultsc');
                for ii=1:inequalitystateconstraintnum.value
                    ocStruct.variable.lagrangemultsc.name{ii}=[basename num2str(ii)];
                end
                ocStruct.variable.lagrangemultsc.num=inequalitystateconstraintnum.value;
            end
        end

        if isfield(ocStruct.objective,'integral')
            if ~isfield(ocStruct.objective.integral,'discountrate') || isempty(ocStruct.objective.integral.discountrate)
                ocStruct.objective.integral.discountrate=getbasicname('discountrate');
            end
            if  ~isfield(ocStruct.objective.integral,'discountfactor') || ~isfield(ocStruct.objective.integral.discountfactor,'type') ...
                    || isempty(ocStruct.objective.integral.discountfactor.type)
                ocStruct.objective.integral.discountfactor.type='expdisc';
            end

            %
            if strcmpi(ocStruct.objective.integral.discountfactor.type,'expdisc')
                ocStruct.objective.integral.discountfactor.term=['exp(-' ocStruct.objective.integral.discountrate '*' ocStruct.variable.independent.name ')'];
            end
        end

        salvagevalue=retrievemodelinformation(ocStruct,'salvagevalue');
        if ~isempty(salvagevalue.value)
            if ~isfield(ocStruct.objective,'sum')
                ocStruct.objective.sum.endtime.function.term=char(salvagevalue.value);
            end
            if ~isfield(ocStruct.objective.sum.endtime,'discountrate') || isempty(ocStruct.objective.sum.endtime.discountrate)
                if isfield(ocStruct.objective.integral,'discountrate') && ~isempty(ocStruct.objective.integral.discountrate)
                    ocStruct.objective.sum.endtime.discountrate=ocStruct.objective.integral.discountrate;
                else
                    ocStruct.objective.sum.endtime.discountrate=getbasicname('discountrate');
                end
            end
            if  ~isfield(ocStruct.objective.sum.endtime,'discountfactor') || ~isfield(ocStruct.objective.sum.endtime.discountfactor,'type') ...
                    || isempty(ocStruct.objective.sum.endtime.discountfactor.type)
                ocStruct.objective.sum.endtime.discountfactor.type='expdisc';
            end

            %
            if strcmpi(ocStruct.objective.sum.endtime.discountfactor.type,'expdisc')
                ocStruct.objective.sum.endtime.discountfactor.term=['exp(-' ocStruct.objective.integral.discountrate '*' ocStruct.variable.endtime.name ')'];
            end
        end

        % if optimization type is not set 'max' is assumed
%         if ~isfield(ocStruct,'optimizationtype')
%             ocStruct.optimizationtype='max';
%         end
        if ~isfield(ocStruct,'optimization')
            ocStruct.optimization.type='max';
            ocStruct.optimization.method='bvp';
        end

        % setting the control value properties 'linear' vs 'nonlinear' and
        % 'explicit' vs 'implicit'
        controlnum=retrievemodelinformation(ocStruct,'controlnum');
        arcnum=retrievemodelinformation(ocStruct,'arcnum');
        arcidentifier=retrievemodelinformation(ocStruct,'arcidentifier');
        defaultexplicitval=num2cell(ones(1,controlnum.value));
        defaultlinearval=num2cell(zeros(1,controlnum.value));
        for ii=1:arcnum.value
            arc=arcidentifier2field(arcidentifier.value{ii});
            try
                propertyStruct=ocStruct.variable.control.(arc).property;
                for propertyname={'explicit','linear'}
                    if ~isfield(propertyStruct,propertyname{1})
                        propertyStruct.(propertyname{1})=eval(['default' propertyname{1} 'val']);
                    end
                end
            catch
                propertyStruct.explicit=defaultexplicitval;
                propertyStruct.linear=defaultlinearval;
            end
            ocStruct.variable.control.(arc).property=propertyStruct;
            clear propertyStruct;
        end

        % if variationparameter exists add variation variables to the list
        % of state and costate variables
        if isfield(ocStruct.variable,'variationparameter')
            ctr1=0;
            ctr2=0;
            ctr3=0;
            for ii=1:ocStruct.variable.variationparameter.num
                for jj=1:ocStruct.variable.state.num
                    ctr1=ctr1+1;
                    ocStruct.variable.variationstate.name{ctr1}=[ocStruct.variable.state.name{jj} ocStruct.variable.variationparameter.name{ii}];
                end
                for jj=1:ocStruct.variable.costate.num
                    ctr2=ctr2+1;
                    ocStruct.variable.variationcostate.name{ctr2}=[ocStruct.variable.costate.name{jj} ocStruct.variable.variationparameter.name{ii}];
                end
                for jj=1:ocStruct.variable.control.num
                    ctr3=ctr3+1;
                    ocStruct.variable.variationcontrol.name{ctr3}=[ocStruct.variable.control.name{jj} ocStruct.variable.variationparameter.name{ii}];
                end
            end
            ocStruct.variable.variationstate.num=length(ocStruct.variable.variationstate.name);
            ocStruct.variable.variationcostate.num=length(ocStruct.variable.variationcostate.name);
            ocStruct.variable.variationcontrol.num=length(ocStruct.variable.variationcontrol.name);
        end

        % determine if the problem is autonomous, neither the objective function
        % nor the state dynamics depend explicitly on the independent variable
        ocStruct.variable.independent.property.autonomous=testautonomous(ocStruct);

    case 'focinitfile'
        arcnum=retrievemodelinformation(ocStruct,'arcnum');
        arcidentifier=retrievemodelinformation(ocStruct,'arcidentifier');
        for ii=1:arcnum.value
            arc=arcidentifier2field(arcidentifier.value{ii});
            for jj=1:ocStruct.variable.control.num
                ocStruct.foc.value.control.(arc).term{jj}=[];
            end
        end

        inequalitycontrolconstraintnum=retrievemodelinformation(ocStruct,'inequalitycontrolconstraintnum');
        if inequalitycontrolconstraintnum.value
            for ii=1:arcnum.value
                arc=arcidentifier2field(arcidentifier.value{ii});
                for jj=1:ocStruct.variable.lagrangemultcc.num
                    ocStruct.foc.value.lagrangemultcc.(arc).term{jj}='0';
                end
            end
        end
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

if isfield(ocStruct.objective,'integral') && ...
        isfield(ocStruct.objective.integral,'discountfactor') && ...
        strcmp(ocStruct.objective.integral.discountfactor.type,'expdisc')
    autflag=isempty(regexp(ocStruct.objective.integral.function.term,['\<' independentvar '\>']));
else
    autflag=1;
end
for ii=1:numel(ocStruct.constraint.derivative.state.term)
    if ~isempty(ocStruct.constraint.derivative.state.term{ii})
        autflag=autflag && isempty(regexp(ocStruct.constraint.derivative.state.term{ii},['\<' independentvar '\>']));
    end
end
if isfield(ocStruct.constraint,'function') && isfield(ocStruct.constraint.function,'control')
    for ii=1:ocStruct.constraint.function.control.num
        if ~isempty(ocStruct.constraint.function.control.term{ii})
            autflag=autflag && isempty(regexp(ocStruct.constraint.function.control.term{ii},['\<' independentvar '\>']));
        end
    end
end
if isfield(ocStruct,'exogenousfunction')
    fn=fieldnames(ocStruct.exogenousfunction);
    for ii=1:length(fn)
            autflag=autflag && isempty(regexp(ocStruct.exogenousfunction.(fn{ii}).term,['\<' independentvar '\>']));
    end
end