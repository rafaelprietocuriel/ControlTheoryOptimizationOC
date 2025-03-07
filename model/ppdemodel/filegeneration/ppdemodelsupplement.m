function ocStruct=ppdemodelsupplement(ocStruct,supplementtype)
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
        
        % set default variable name for independent variables time and
        % space
        if ~isfield(ocStruct.variable,'independent') || ~isfield(ocStruct.variable.independent,'time')
            ocStruct.variable.independent.time.name=getbasicname('time');
        end
        if ~isfield(ocStruct.variable.independent,'space')
            ocStruct.variable.independent.space.name=getbasicname('space');
            ocStruct.variable.independent.space.num=1;
        end
        
        % get time and space variable names
        timearg=retrieveppdemodelinformation(ocStruct,'time');
        spacearg=retrieveppdemodelinformation(ocStruct,'space');
        
        % set the default dependence on the independent variables for
        % state(s) and control(s)
        statenum=retrieveppdemodelinformation(ocStruct,'statenum');
        for ii=1:statenum.value
            if isempty(ocStruct.variable.state(ii).dependence)
                ocStruct.variable.state(ii).dependence='2'; % 0 only time dependent, 1 only space dependent, 2 time and space dependent
            else
                if any(strcmp(ocStruct.variable.state(ii).dependence,timearg.value))
                    if any(strcmp(ocStruct.variable.state(ii).dependence,spacearg.value))
                        ocStruct.variable.state(ii).dependence='2';
                    else
                        ocStruct.variable.state(ii).dependence='0';
                    end
                else
                    ocStruct.variable.state(ii).dependence='1';
                end
            end
            
        end
        
        controlnum=retrieveppdemodelinformation(ocStruct,'controlnum');
        for ii=1:controlnum.value
            if isempty(ocStruct.variable.control(ii).dependence)
                ocStruct.variable.control(ii).dependence='2'; % 0 only time dependent, 1 only space dependent, 2 time and space dependent
            else
                if any(strcmp(ocStruct.variable.control(ii).dependence,timearg.value))
                    if any(strcmp(ocStruct.variable.control(ii).dependence,spacearg.value))
                        ocStruct.variable.control(ii).dependence='2';
                    else
                        ocStruct.variable.control(ii).dependence='0';
                    end
                else
                    ocStruct.variable.control(ii).dependence='1';
                end
            end
            
        end
        
        % setting the control value properties 'linear' vs 'nonlinear' and
        % 'explicit' vs 'implicit'
        arcnum=retrieveppdemodelinformation(ocStruct,'arcnum');
        arcidentifier=retrieveppdemodelinformation(ocStruct,'arcidentifier');
        defaultexplicitval=1;
        defaultlinearval=0;
        for ii=1:arcnum.value
            arc=arcidentifier2field(arcidentifier.value{ii});
            for jj=1:controlnum.value
                try
                    propertyStruct=ocStruct.variable.control(jj).(arc).property;
                    for propertyname={'explicit','linear'}
                        if ~isfield(propertyStruct,propertyname{1})
                            propertyStruct.(propertyname{1})=eval(['default' propertyname{1} 'val']);
                        end
                    end
                catch
                    propertyStruct.explicit=defaultexplicitval;
                    propertyStruct.linear=defaultlinearval;
                end
                ocStruct.variable.control(jj).(arc).property=propertyStruct;
                clear propertyStruct;
            end
        end

        % set default costate name if not set by the user
        if ~isfield(ocStruct.variable,'costate')
            basename=getbasicname('costate');
            for ii=1:statenum.value
                if ~isfield(ocStruct.variable,'costate') ||~isfield(ocStruct.variable.costate(ii),'name') || isempty(ocStruct.variable.costate(ii).name)
                    ocStruct.variable.costate(ii).name=[basename num2str(ii)];
                end
            end
        end
        for ii=1:statenum.value
            ocStruct.variable.costate(ii).dependence=ocStruct.variable.state(ii).dependence;
        end
        
        % set default variable name for end time variable
        if ~isfield(ocStruct.variable,'endtime') || ~isfield(ocStruct.variable.endtime,'name') ...
                || isempty(ocStruct.variable.endtime.name)
            ocStruct.variable.endtime.name=getbasicname('endtime');
            ocStruct.variable.endtime.num=1;
        end

        % set default variable name for Lagrange multiplier if control constraints
        % are given
        inequalitycontrolconstraintnum=retrieveppdemodelinformation(ocStruct,'inequalitycontrolconstraintnum');
        if inequalitycontrolconstraintnum.value
            if ~isfield(ocStruct.variable,'lagrangemultcc') || ~isfield(ocStruct.variable.lagrangemultcc,'name') ...
                    || isempty(ocStruct.variable.lagrangemultcc.name)
                basename=getbasicname('lagrangemultcc');
                for ii=1:inequalitycontrolconstraintnum.value
                    ocStruct.variable.lagrangemultcc(ii).name=[basename num2str(ii)];
                end
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
            discountrate=retrieveppdemodelinformation(ocStruct,'discountrate');
            
            if strcmpi(ocStruct.objective.integral.discountfactor.type,'expdisc')
                ocStruct.objective.integral.discountfactor.term=['exp(-' discountrate.value '*' timearg.value ')'];
            end
            
            if ~isfield(ocStruct.objective.integral,'dependence') || isempty(ocStruct.objective.integral.dependence)
                ocStruct.objective.integral.function.dependence='2';
            else
                if any(strcmp(ocStruct.objective.integral.dependence,timearg.value))
                    if any(strcmp(ocStruct.objective.integral.function.dependence,spacearg.value))
                        ocStruct.objective.integral.function.dependence='2';
                    else
                        ocStruct.objective.integral.function.dependence='0';
                    end
                else
                    ocStruct.objective.integral.function.dependence='1';
                end
            end
        end

        salvagevalue=retrieveppdemodelinformation(ocStruct,'salvagevalue');
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
        if ~isfield(ocStruct,'optimizationtype')
            ocStruct.optimizationtype='max';
        end

        % determine if the problem is autonomous, neither the objective function
        % nor the state dynamics depend explicitly on the independent variable
        ocStruct.variable.independent.property.autonomous=testautonomous(ocStruct);

        % determine if the problem explicitly depends on the spatial variable, neither the objective function
        % nor the state dynamics depend explicitly on the independent variable
        ocStruct.variable.independent.property.explicitspatial=testexplicitspatial(ocStruct);

        femdataXcoord=retrieveppdemodelinformation(ocStruct,'femdataXcoord');

        femdataXcoordm1=retrieveppdemodelinformation(ocStruct,'femdataXcoordm1');

        femdataXcoordp1=retrieveppdemodelinformation(ocStruct,'femdataXcoordp1');
        
        femdataucoord=retrieveppdemodelinformation(ocStruct,'femdataucoord');

        ocStruct.spacegeometry.femdata.Xcoord=femdataXcoord.value;
        ocStruct.spacegeometry.femdata.Xcoordm1=femdataXcoordm1.value;
        ocStruct.spacegeometry.femdata.Xcoordp1=femdataXcoordp1.value;
        ocStruct.spacegeometry.femdata.ucoord=femdataucoord.value;

    case 'focinitfile'
        arcnum=retrieveppdemodelinformation(ocStruct,'arcnum');
        arcidentifier=retrieveppdemodelinformation(ocStruct,'arcidentifier');
        for ii=1:arcnum.value
            arc=arcidentifier2field(arcidentifier.value{ii});
            for jj=1:ocStruct.variable.control.num
                ocStruct.foc.value.control.(arc).term{jj}=[];
            end
        end

        inequalitycontrolconstraintnum=retrieveppdemodelinformation(ocStruct,'inequalitycontrolconstraintnum');
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
timearg=retrieveppdemodelinformation(ocStruct,'time');
inequalitycontrolconstraintnum=retrieveppdemodelinformation(ocStruct,'inequalitycontrolconstraintnum');

if isfield(ocStruct.objective,'integral') && ...
        isfield(ocStruct.objective.integral,'discountfactor') && ...
        strcmp(ocStruct.objective.integral.discountfactor.type,'expdisc')
    autflag=isempty(regexp(ocStruct.objective.integral.function.term,['\<' timearg.value '\>']));
else
    autflag=1;
end
for ii=1:numel(ocStruct.constraint.derivative.state)
    if ~isempty(ocStruct.constraint.derivative.state(ii).term)
        autflag=autflag && isempty(regexp(ocStruct.constraint.derivative.state(ii).term,['\<' timearg.value '\>']));
    end
end
if isfield(ocStruct.constraint,'function') && isfield(ocStruct.constraint.function,'control')
    for ii=1:inequalitycontrolconstraintnum.value
        if ~isempty(ocStruct.constraint.function.control(ii).term)
            autflag=autflag && isempty(regexp(ocStruct.constraint.function.control(ii).term,['\<' timearg.value '\>']));
        end
    end
end
if isfield(ocStruct,'exogenousfunction')
    fn=fieldnames(ocStruct.exogenousfunction);
    for ii=1:length(fn)
            autflag=autflag && isempty(regexp(ocStruct.exogenousfunction.(fn{ii}).term,['\<' timearg.value '\>']));
    end
end


function spatialflag=testexplicitspatial(ocStruct)
%
% the OC problem is called autonomous if the integrand of the objective
% function, state dynamics and constraint functions do not explicitly
% depend on time. Furthermore it is assumed that the discount factor is
% given by the exponential function. Under these assumptions the canonical
% system can be formulated as an autonomous system of ODEs justifying the
% terminology
spatialflag=0;
if isempty(ocStruct)
    return
end
spacearg=retrieveppdemodelinformation(ocStruct,'space');
inequalitycontrolconstraintnum=retrieveppdemodelinformation(ocStruct,'inequalitycontrolconstraintnum');

if isfield(ocStruct.objective,'integral') && ...
        isfield(ocStruct.objective.integral,'discountfactor') && ...
        strcmp(ocStruct.objective.integral.discountfactor.type,'expdisc')
    spatialflag=isempty(regexp(ocStruct.objective.integral.function.term,['\<' spacearg.value '\>']));
else
    spatialflag=1;
end
for ii=1:numel(ocStruct.constraint.derivative.state)
    if ~isempty(ocStruct.constraint.derivative.state(ii).term)
        spatialflag=spatialflag && isempty(regexp(ocStruct.constraint.derivative.state(ii).term,['\<' spacearg.value '\>']));
    end
end
if isfield(ocStruct.constraint,'function') && isfield(ocStruct.constraint.function,'control')
    for ii=1:inequalitycontrolconstraintnum.value
        if ~isempty(ocStruct.constraint.function.control(ii).term)
            spatialflag=spatialflag && isempty(regexp(ocStruct.constraint.function.control(ii).term,['\<' spacearg.value '\>']));
        end
    end
end
if isfield(ocStruct,'exogenousfunction')
    fn=fieldnames(ocStruct.exogenousfunction);
    for ii=1:length(fn)
            spatialflag=spatialflag && isempty(regexp(ocStruct.exogenousfunction.(fn{ii}).term,['\<' spacearg.value '\>']));
    end
end
spatialflag=~spatialflag;