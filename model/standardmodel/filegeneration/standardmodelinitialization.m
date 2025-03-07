function ocStruct=standardmodelinitialization(initsection,initfiletext,ocStruct)
%
%
propname=fieldnames(initsection);
for ii=1:numel(propname)
    try
        ocStruct=propertyhandling(ocStruct,initfiletext(initsection.(propname{ii})(1):initsection.(propname{ii})(2)),propname{ii});
    catch
        ocStruct=[];
        ocmatmsg('Initialization process aborted.\nProblems in reading section ''%s''.\n',propname{ii})
        return
    end
end

function ocStruct=propertyhandling(ocStruct,sectiontext,propname)

switch propname
    case 'description'
        ocStruct.description=sectiontext;
    case 'variable'
        for ii=1:numel(sectiontext)
            [leftstr,rightstr]=getpropstrings(sectiontext{ii});
            leftstr=[leftstr{:}];
            switch leftstr
                case 'independent'
                    ocStruct.variable.independent.name=[rightstr{:}];
                case 'endtime'
                    ocStruct.variable.(leftstr).name=[rightstr{:}];
                case 'variationparameter'
                    ocStruct.variable.variationparameter.name=rightstr;
                    ocStruct.variable.variationparameter.num=numel(rightstr);
                otherwise
                    ocStruct.variable.(leftstr).name=rightstr;
                    if strcmp(rightstr,'[]')
                        ocStruct.variable.(leftstr).num=0;
                    else
                        ocStruct.variable.(leftstr).num=numel(rightstr);
                    end
            end
        end
    case 'control'
        for ii=1:numel(sectiontext)
            [leftstr,middlestr,rightstr]=getpropstrings(sectiontext{ii});
            if strcmp([leftstr{:}],'*')
                leftstr=ocStruct.arc.identifier;
            end
            if strcmp([middlestr{:}],'*')
                middlestr=ocStruct.variable.control.name;
            end
            if ~any(strcmp(rightstr,'replace'))
                rightstr=[rightstr{:}];
                switch rightstr
                    case 'explicit'
                        for jj=1:numel(leftstr) % arc identifiers
                            for kk=1:numel(middlestr)
                                ctrlidx=find(strcmp(middlestr{kk},ocStruct.variable.control.name));
                                ocStruct.variable.control.(arcidentifier2field(leftstr{jj})).property.explicit{ctrlidx}=1;
                            end
                        end
                    case 'implicit'
                        for jj=1:numel(leftstr) % arc identifiers
                            for kk=1:numel(middlestr)
                                ctrlidx=find(strcmp(middlestr{kk},ocStruct.variable.control.name));
                                ocStruct.variable.control.(arcidentifier2field(leftstr{jj})).property.explicit{ctrlidx}=0;
                            end
                        end
                    case 'nonlinear'
                        for jj=1:numel(leftstr) % arc identifiers
                            for kk=1:numel(middlestr)
                                ctrlidx=find(strcmp(middlestr{kk},ocStruct.variable.control.name));
                                ocStruct.variable.control.(arcidentifier2field(leftstr{jj})).property.linear{ctrlidx}=0;
                            end
                        end
                    case 'linear'
                        for jj=1:numel(leftstr) % arc identifiers
                            for kk=1:numel(middlestr)
                                ctrlidx=find(strcmp(middlestr{kk},ocStruct.variable.control.name));
                                ocStruct.variable.control.(arcidentifier2field(leftstr{jj})).property.linear{ctrlidx}=1;
                            end
                        end
                end
            else
                rightstr(strcmp(rightstr,'replace'))=[];
                for jj=1:numel(leftstr) % arc identifiers
                    for kk=1:numel(middlestr)
                        ctrlidx=find(strcmp(middlestr{kk},ocStruct.variable.control.name));
                        ocStruct.variable.control.(arcidentifier2field(leftstr{jj})).property.replace{ctrlidx}=str2num(rightstr{1});
                    end
                end
            end
        end
        
    case 'state'
        for ii=1:numel(sectiontext)
            [leftstr,middlestr]=getpropstrings(sectiontext{ii});
            for jj=1:numel(leftstr) % arc identifiers
                ocStruct.variable.state.(arcidentifier2field(leftstr{jj})).explicit.name=middlestr;
            end
        end
        

    case 'statedynamics'
        odecounter=0;
        statevarname=ocStruct.variable.state.name;
        for ii=1:numel(sectiontext)
            [leftstr,rightstr]=getpropstrings(sectiontext{ii});
            leftstr=[leftstr{:}];
            rightstr=[rightstr{:}];
            switch leftstr
                case 'ode'
                    odecounter=odecounter+1;
                    fidx=strfind(rightstr,'=');
                    leftidx=[1 fidx+1];
                    rightidx=[fidx-1 numel(rightstr)];
                    derivative=strtrim(rightstr(leftidx(1):rightidx(1)));
                    fidx=find(strcmp(statevarname,derivative(2:end)));
                    ocStruct.constraint.derivative.state.type{fidx}=leftstr;
                    ocStruct.constraint.derivative.state.term{fidx}=strtrim(rightstr(leftidx(2):rightidx(2)));
            end
        end
        ocStruct.constraint.derivative.state.num=odecounter;

    case 'objective'
        for ii=1:numel(sectiontext)
            [leftstr,rightstr]=getpropstrings(sectiontext{ii});
            leftstr=[leftstr{:}];
            rightstr=[rightstr{:}];
            switch leftstr
                case 'int'
                    ocStruct.objective.integral.function.term=rightstr;
                case 'expdisc'
                    ocStruct.objective.integral.discountrate=rightstr;
                    ocStruct.objective.integral.discountfactor.type=leftstr;
            end
        end
        
    case 'salvagevalue'
        for ii=1:numel(sectiontext)
            [leftstr,rightstr]=getpropstrings(sectiontext{ii});
            leftstr=[leftstr{:}];
            rightstr=[rightstr{:}];
            switch leftstr
                case 'func'
                    ocStruct.objective.sum.endtime.function.term=rightstr;
                case 'expdisc'
                    ocStruct.objective.sum.endtime.discountrate=rightstr;
                    %ocStruct.objective.sum.endtime.discountfactor.term=['exp(-' rightstr '*' ocStruct.variable.endtime.name ')'];
            end
        end
        
    case 'optimization'
        for ii=1:numel(sectiontext)
            [leftstr,rightstr]=getpropstrings(sectiontext{ii});
            leftstr=[leftstr{:}];
            switch leftstr
                case 'type'
                    ocStruct.optimization.type.term=[rightstr{:}];
                case 'method'
                    ocStruct.optimization.method=char(rightstr);
            end
        end
        
    case 'optimizationtype'
        ocStruct.optimization.type=sectiontext{1};

    case 'controlconstraint'
        for ii=1:numel(sectiontext)
            [leftstr,middlestr,rightstr]=getpropstrings(sectiontext{ii});
            middlestr=[middlestr{:}];
            leftstr=[leftstr{:}];
            rightstr=[rightstr{:}];
            ocStruct.constraint.function.control.identifier{ii}=leftstr;
            ocStruct.constraint.function.control.type{ii}=middlestr;
            switch middlestr
                case 'ineq'
                    fidx=findstr(rightstr,'>=');
                    if ~isempty(fidx)
                        leftidx=[1 fidx+2];
                        rightidx=[fidx-1 numel(rightstr)];
                        ocStruct.constraint.function.control.term{ii}=[strtrim(rightstr(leftidx(1):rightidx(1))) '-(' strtrim(rightstr(leftidx(2):rightidx(2))) ')'];
                    end
                    fidx=findstr(rightstr,'<=');
                    if ~isempty(fidx)
                        leftidx=[1 fidx+2];
                        rightidx=[fidx-1 numel(rightstr)];
                        ocStruct.constraint.function.control.term{ii}=[strtrim(rightstr(leftidx(2):rightidx(2))) '-(' strtrim(rightstr(leftidx(1):rightidx(1))) ')'];
                    end
            end
        end
        ocStruct.constraint.function.control.num=ii;
        

    case 'stateconstraint'
        for ii=1:numel(sectiontext)
            [leftstr,middlestr,rightstr]=getpropstrings(sectiontext{ii});
            middlestr=[middlestr{:}];
            leftstr=[leftstr{:}];
            ocStruct.constraint.function.state.identifier{ii}=leftstr;
            ocStruct.constraint.function.state.type{ii}=middlestr;
            switch middlestr
                case 'ineq'
                    lenrs=length(rightstr);
                    if lenrs==2
                        if strfind(rightstr{1},'=')
                            order=rightstr{2};
                            rightstr=rightstr{1};
                        else
                            order=rightstr{1};
                            rightstr=rightstr{2};
                        end
                    else
                        ocmaterror('Order of state constraint has to be provided.')
                    end
                    ocStruct.constraint.function.state.order{ii}=str2double(order);
                    fidx=findstr(rightstr,'>=');
                    if ~isempty(fidx)
                        leftidx=[1 fidx+2];
                        rightidx=[fidx-1 numel(rightstr)];
                        ocStruct.constraint.function.state.term{ii}=[strtrim(rightstr(leftidx(1):rightidx(1))) '-(' strtrim(rightstr(leftidx(2):rightidx(2))) ')'];
                    end
                    fidx=findstr(rightstr,'<=');
                    if ~isempty(fidx)
                        leftidx=[1 fidx+2];
                        rightidx=[fidx-1 numel(rightstr)];
                        ocStruct.constraint.function.state.term{ii}=[strtrim(rightstr(leftidx(2):rightidx(2))) '-(' strtrim(rightstr(leftidx(1):rightidx(1))) ')'];
                    end
            end
        end
        ocStruct.constraint.function.state.num=ii;
        
    case 'arcdefinition'
        for ii=1:numel(sectiontext)
            [leftstr,rightstr]=getpropstrings(sectiontext{ii});
            leftstr=[leftstr{:}];
            %rightstr=	;
            rightstr=sprintf('%s ',rightstr{:});
            rightstr(end)=[];
            ocStruct.arc.identifier{ii}=leftstr;
            ocStruct.arc.argument(ii)=str2double(leftstr);
            fidx=findstr(rightstr,'_');
            if ~isempty(fidx)
                leftidx=[1 fidx+1];
                rightidx=[fidx-1 numel(rightstr)];
                ocStruct.arc.constraintcombination{ii}=rightstr(leftidx(1):rightidx(1));
                ocStruct.arc.controlvaluecombination{ii}=rightstr(leftidx(2):rightidx(2));
            else
                ocStruct.arc.constraintcombination{ii}=rightstr;
                ocStruct.arc.controlvaluecombination{ii}='1';
            end

        end
        ocStruct.arc.num=ii;
        
        
    case {'maximizingderivativevariable','maximizingexplicitvariable'}
        for ii=1:numel(sectiontext)
            [leftstr,rightstr]=getpropstrings(sectiontext{ii});
            leftstr=[leftstr{:}];
            ocStruct.foc.generalinformation.(arcidentifier2field(leftstr)).(propname).name=rightstr;
        end
        
    case 'parameter'
        ocStruct.parameter.variable=[];
        ocStruct.parameter.algebraictermidx=[];
        countparameter=0;
        for ii=1:numel(sectiontext)
            [leftstr,rightstr]=getpropstrings(sectiontext{ii});
            leftstr=[leftstr{:}];
            rightstr=[rightstr{:}];
            try
                if ~isfield(ocStruct.parameter.variable,leftstr)
                    ocStruct.parameter.variable.(leftstr)=eval(rightstr);
                    countparameter=countparameter+1;
                end
            catch
                if ~isfield(ocStruct.parameter.variable,leftstr)
                    countparameter=countparameter+1;
                    ocStruct.parameter.variable.(leftstr)=rightstr;
                    ocStruct.parameter.algebraictermidx=[ocStruct.parameter.algebraictermidx ii];
                end
            end
        end
        ocStruct.parameter.num=countparameter;

    case 'exogenousfunction'
        for ii=1:numel(sectiontext)
            [leftstr,rightstr]=getpropstrings(sectiontext{ii});
            leftstr=[leftstr{:}];
            rightstr=[rightstr{:}];
            argidx=regexp(leftstr, '\(','once');
            if ~isempty(argidx)
                name=leftstr(1:argidx-1);
                ocStruct.exogenousfunction.(name).argument=leftstr(argidx+1:end-1);
            else
                name=leftstr;
                ocStruct.exogenousfunction.(name).argument='';
            end
            
            if ~isempty(rightstr)
                ocStruct.exogenousfunction.(name).term=rightstr;
            else
                ocStruct.exogenousfunction.(name).term='';
            end
        end

    case 'exogenousdynamics'
        statevarname=ocStruct.variable.exogenousstate.name;
        odecounter=0;
        for ii=1:numel(sectiontext)
            [leftstr,rightstr]=getpropstrings(sectiontext{ii});
            leftstr=[leftstr{:}];
            rightstr=[rightstr{:}];
            switch leftstr
                case 'ode'
                    odecounter=odecounter+1;
                    fidx=strfind(rightstr,'=');
                    leftidx=[1 fidx+1];
                    rightidx=[fidx-1 numel(rightstr)];
                    derivative=strtrim(rightstr(leftidx(1):rightidx(1)));
                    fidx=find(strcmp(statevarname,derivative(2:end)));
                    ocStruct.exogenousdynamics.derivative.state.type{fidx}=leftstr;
                    ocStruct.exogenousdynamics.derivative.state.term{fidx}=strtrim(rightstr(leftidx(2):rightidx(2)));
            end
        end
%         ocStruct.constraint.derivative.state.num=odecounter;
%         for ii=1:numel(sectiontext)
%             [leftstr,rightstr]=getpropstrings(sectiontext{ii});
%             leftstr=[leftstr{:}];
%             rightstr=[rightstr{:}];
%             argidx=regexp(leftstr, '\(','once');
%             if ~isempty(argidx)
%                 name=leftstr(1:argidx-1);
%                 ocStruct.exogenousdynamics.(name).argument=leftstr(argidx+1:end-1);
%             else
%                 name=leftstr;
%                 ocStruct.exogenousdynamics.(name).argument='';
%             end
%             
%             if ~isempty(rightstr)
%                 ocStruct.exogenousdynamics.(name).term=rightstr;
%             else
%                 ocStruct.exogenousdynamics.(name).term='';
%             end
%         end

    case 'nonsmoothfunction'
        for ii=1:numel(sectiontext)
            [leftstr,middlestr,rightstr]=getpropstrings(sectiontext{ii});
            middlestr=[middlestr{:}];
            leftstr=[leftstr{:}];
            rightstr=[rightstr{:}];
            argidx=regexp(leftstr, '\(','once');
            if ~isempty(argidx)
                name=leftstr(1:argidx-1);
            else
                name=leftstr;
            end
            identifier=middlestr;
            if ~isempty(argidx)
                ocStruct.nonsmoothfunction.(name).(identifier).argument=leftstr(argidx+1:end-1);
            else
                ocStruct.nonsmoothfunction.(name).(identifier).argument='';
            end

            if ~isempty(rightstr)
                ocStruct.nonsmoothfunction.(name).(identifier).term=rightstr;
            else
                ocStruct.nonsmoothfunction.(name).(identifier).term='';
            end
            ocStruct.nonsmoothfunction.(name).(identifier).name=name;
        end

    case 'switchingcondition'
        for ii=1:numel(sectiontext)
            [leftstr,middlestr,rightstr]=getpropstrings(sectiontext{ii});
            middlestr=[middlestr{:}];
            leftstr=[leftstr{:}];
            rightstr=[rightstr{:}];
            ocStruct.switchingcondition.identifier{ii}=leftstr;
            ocStruct.switchingcondition.type{ii}=middlestr;
            switch middlestr
                case 'ineq'
                    fidx=findstr(rightstr,'>=');
                    if ~isempty(fidx)
                        leftidx=[1 fidx+2];
                        rightidx=[fidx-1 numel(rightstr)];
                        ocStruct.switchingcondition.term{ii}=[strtrim(rightstr(leftidx(1):rightidx(1))) '-(' strtrim(rightstr(leftidx(2):rightidx(2))) ')'];
                    end
                    fidx=findstr(rightstr,'<=');
                    if ~isempty(fidx)
                        leftidx=[1 fidx+2];
                        rightidx=[fidx-1 numel(rightstr)];
                        ocStruct.switchingcondition.term{ii}=[strtrim(rightstr(leftidx(2):rightidx(2))) '-(' strtrim(rightstr(leftidx(1):rightidx(1))) ')'];
                    end
            end
        end
        ocStruct.switchingcondition.num=ii;
end
function varargout=getpropstrings(inputstr)
inputstr=regexprep(inputstr,'\s',''); % remove all white-space characters
varargout=cell(nargout,1);
tmpcell=regexp(inputstr,'::','split');
if ~isempty(tmpcell)
    varargout(1:numel(tmpcell))=tmpcell;
else
    varargout{1}='';
    varargout{2}=inputstr;
end
% if output argument is comma or colon separated transform it into cell
% array
for ii=1:numel(varargout)
    if ~isempty(strfind(varargout{ii},','))
        if prod(size(mystr2sym(varargout{ii})))>1
            varargout{ii}=regexp(varargout{ii},',','split');
        else
            varargout{ii}={varargout{ii}};
        end
    elseif numel(strfind(varargout{ii},':'))==1
        tmpcell=regexp(varargout{ii},':','split');
        varname=regexp(tmpcell,'[a-zA-Z_]*','match');
        if ~strcmp(varname{1},varname{2})
            ocmaterror('Initialization aborted. Variable names have to be the same with a colon operator.')
        end
        varidx=regexp(tmpcell,'[0-9]*','match');
        idx=[str2double(varidx{1}):str2double(varidx{2})];
        varargout{ii}=cell(1,length(idx));
        ctr=0;
        for jj=idx
            ctr=ctr+1;
            if ~isempty(varname{1})
                varargout{ii}{ctr}=[char(varname{1}) num2str(jj)];
            else
                varargout{ii}{ctr}=[num2str(jj)];
            end
        end
    else
        varargout{ii}={varargout{ii}};
    end
end
