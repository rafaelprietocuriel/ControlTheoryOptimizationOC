function ocStruct=standarddiffmodelinitialization(initsection,initfiletext,ocStruct)
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
            [leftstr rightstr]=getpropstrings(sectiontext{ii});
            leftstr=[leftstr{:}];
            switch leftstr
                case 'independent'
                    ocStruct.variable.independent.name=[rightstr{:}];
                case 'endtime'
                    ocStruct.variable.(leftstr).name=[rightstr{:}];
                otherwise
                    ocStruct.variable.(leftstr).name=rightstr;
                    ocStruct.variable.(leftstr).num=numel(rightstr);
            end
        end
    case 'control'
        for ii=1:numel(sectiontext)
            [leftstr middlestr rightstr]=getpropstrings(sectiontext{ii});
            rightstr=[rightstr{:}];
            if strcmp([leftstr{:}],'*')
                leftstr=ocStruct.arc.identifier;
            end
            if strcmp([middlestr{:}],'*')
                middlestr=ocStruct.variable.control.name;
            end
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
        end
        
    case 'statedynamics'
        diffcounter=0;
        statevarname=ocStruct.variable.state.name;
        for ii=1:numel(sectiontext)
            [leftstr middlestr rightstr]=getpropstrings(sectiontext{ii});
            leftstr=[leftstr{:}];
            middlestr=[middlestr{:}];
            rightstr=[rightstr{:}];
            if isempty(rightstr)
                rightstr=middlestr;
                switch leftstr
                    case 'diff'
                        diffcounter=diffcounter+1;
                        fidx=strfind(rightstr,'=');
                        leftidx=[1 fidx+1];
                        rightidx=[fidx-1 numel(rightstr)];
                        statetp1=regexprep(strtrim(rightstr(leftidx(1):rightidx(1))),'_(\w)*','');
                        fidx=find(strcmp(statevarname,statetp1));
                        ocStruct.constraint.difference.state.type{fidx}=[leftstr '_explicit'];
                        ocStruct.constraint.difference.state.term{fidx}=strtrim(rightstr(leftidx(2):rightidx(2)));
                end
            else
            end
        end
        ocStruct.constraint.difference.state.num=diffcounter;

    case 'objective'
        for ii=1:numel(sectiontext)
            [leftstr rightstr]=getpropstrings(sectiontext{ii});
            leftstr=[leftstr{:}];
            rightstr=[rightstr{:}];
            switch leftstr
                case 'sum'
                    ocStruct.objective.sum.function.term=rightstr;
                case 'expdisc'
                    ocStruct.objective.sum.discountrate=rightstr;
                    ocStruct.objective.sum.discountfactor.type=leftstr;
            end
        end
        
    case 'salvagevalue'
        for ii=1:numel(sectiontext)
            [leftstr rightstr]=getpropstrings(sectiontext{ii});
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
        
    case 'optimizationtype'
        ocStruct.optimizationtype=sectiontext;

    case 'controlconstraint'
        for ii=1:numel(sectiontext)
            [leftstr middlestr rightstr]=getpropstrings(sectiontext{ii});
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
        
    case 'arcdefinition'
        for ii=1:numel(sectiontext)
            [leftstr rightstr]=getpropstrings(sectiontext{ii});
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
        
        
    case 'maptype'
        ocStruct.constraint.difference.maptype=sectiontext;

    case {'maximizingderivativevariable','maximizingexplicitvariable'}
        for ii=1:numel(sectiontext)
            [leftstr rightstr]=getpropstrings(sectiontext{ii});
            leftstr=[leftstr{:}];
            ocStruct.foc.generalinformation.(arcidentifier2field(leftstr)).(propname).name=rightstr;
        end
        
    case 'parameter'
        for ii=1:numel(sectiontext)
            [leftstr rightstr]=getpropstrings(sectiontext{ii});
            leftstr=[leftstr{:}];
            rightstr=[rightstr{:}];
            ocStruct.parameter.variable.(leftstr)=eval(rightstr);
        end
        ocStruct.parameter.num=ii;
        
    case 'exogenousfunction'
        for ii=1:numel(sectiontext)
            [leftstr rightstr]=getpropstrings(sectiontext{ii});
            leftstr=[leftstr{:}];
            rightstr=[rightstr{:}];
            ocStruct.exogenousfunction.(leftstr).term=rightstr;
        end
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
        varargout{ii}=regexp(varargout{ii},',','split');
    elseif numel(strfind(varargout{ii},':'))==1
        tmpcell=regexp(varargout{ii},':','split');
        varname=regexp(tmpcell,'[a-zA-Z_]*','match');
        if ~strcmp(varname{1},varname{2})
            ocmaterror('Initialization aborted. Variable names have to be the same with a colon operator.')
        end
        varidx=regexp(tmpcell,'[0-9]*','match');
        idx=[str2double(varidx{1}):str2double(varidx{2})];
        varargout{ii}=cell(1,length(idx));
        for jj=idx
            varargout{ii}{jj}=[char(varname{1}) num2str(jj)];
        end
    else
        varargout{ii}={varargout{ii}};
    end
end
