function ocStruct=ppdemodelinitialization(initsection,initfiletext,ocStruct)
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
            [leftstr,midstr,rightstr]=getpropstrings(sectiontext{ii});
            leftstr=[leftstr{:}];
            if isempty([rightstr{:}])
                rightstr=midstr;
                midstr=[];
            end
            switch leftstr
                case 'time'
                    ocStruct.variable.independent.time.name=[rightstr{:}];
                case 'space'
                    if length(rightstr)==1
                        ocStruct.variable.independent.space.name=[rightstr{:}];
                    else
                        for jj=1:length(rightstr)
                            ocStruct.variable.independent.space(ii).name=rightstr{ii};
                        end
                    end
                case 'endtime'
                    ocStruct.variable.(leftstr).name=[rightstr{:}];
                otherwise
                    if ~isfield(ocStruct,'variable') || ...
                            ~isfield(ocStruct.variable,leftstr) || ...
                            ~isfield(ocStruct.variable.(leftstr),'name')
                        if length(rightstr)==1
                            ocStruct.variable.(leftstr).name=[rightstr{:}];
                            ocStruct.variable.(leftstr).dependence=midstr;
                        else
                            for jj=1:length(rightstr)
                                ocStruct.variable.(leftstr)(jj).name=rightstr{jj};
                                ocStruct.variable.(leftstr)(jj).dependence=midstr;
                            end
                        end
                    else
                        offset=length(ocStruct.variable.(leftstr));
                        for jj=1:length(rightstr)
                            ocStruct.variable.(leftstr)(jj+offset).name=rightstr{jj};
                            ocStruct.variable.(leftstr)(jj+offset).dependence=midstr;
                        end
                    end
            end
        end

    case 'statedynamics'
        odecounter=0;
        pdecounter=0;
        statename=retrieveppdemodelinformation(ocStruct,'statename');
        for ii=1:numel(sectiontext)
            [leftstr,rightstr]=getpropstrings(sectiontext{ii});
            leftstr=[leftstr{:}];
            rightstr=[rightstr{:}];
            switch leftstr
                case {'ode','pde'}
                    if strcmp(leftstr,'pde')
                        pdecounter=pdecounter+1;
                    else
                        odecounter=odecounter+1;
                    end
                    fidx=strfind(rightstr,'=');
                    nonlinearterm=rightstr(fidx+1:end);
                    % find time variable after derivative operator Dt.
                    timederivativevariable=regexpi(rightstr,'(?<=Dt.)\w+','match');
                    % find first order space derivative operator D1.
                    [convectioncoefficient,idx0,dum]=regexpi(nonlinearterm,'[+-]{1}\w+\*(?=D1\.)|[+-]{1}\[\S+\]\*(?=D1\.)','match');
                    [convectionvariable,dum,idx1]=regexpi(nonlinearterm,'(?<=D1\.)\w+|(?<=D1\.)(\[){1}\S*(\]){1}','match');
                    nonlinearterm(idx0:idx1)=[];
                    % find second order space derivative operator D2.
                    [diffusioncoefficient,idx0,dum]=regexpi(nonlinearterm,'[+-]{1}\w+\*(?=D2\.)|[+-]{1}\[\S+\]\*(?=D2\.)','match');
                    [diffusionvariable,dum,idx1]=regexpi(nonlinearterm,'(?<=D2\.)\w+|(?<=D2\.)(\[){1}\S*(\]){1}','match');
                    nonlinearterm(idx0:idx1)=[];

                    fidx=find(strcmp(statename.value,timederivativevariable));
                    ocStruct.constraint.derivative.state(fidx).type=leftstr;
                    ocStruct.constraint.derivative.state(fidx).term=strtrim(nonlinearterm);
                    if ~isempty(convectioncoefficient)
                        convectioncoefficient=[convectioncoefficient{:}];
                        if strcmp(convectioncoefficient(end),'*')
                            convectioncoefficient(end)=[];
                        end
                        if strcmp(convectioncoefficient(1),'+')
                            convectioncoefficient(1)=[];
                        end
                        ocStruct.constraint.derivative.state(fidx).convection.term=convectioncoefficient;
                    end
                    if ~isempty(convectionvariable)
                        ocStruct.constraint.derivative.state(fidx).convection.variable=[convectionvariable{:}];
                    end
                    if ~isempty(diffusioncoefficient)
                        diffusioncoefficient=[diffusioncoefficient{:}];
                        if strcmp(diffusioncoefficient(end),'*')
                            diffusioncoefficient(end)=[];
                        end
                        if strcmp(diffusioncoefficient(1),'+')
                            diffusioncoefficient(1)=[];
                        end
                        ocStruct.constraint.derivative.state(fidx).diffusion.term=diffusioncoefficient;
                    end
                    if ~isempty(diffusionvariable)
                        ocStruct.constraint.derivative.state(fidx).diffusion.variable=[diffusionvariable{:}];
                    end
            end
        end

    case 'objective'
        for ii=1:numel(sectiontext)
            [leftstr,midstr,rightstr]=getpropstrings(sectiontext{ii});
            if isempty([rightstr{:}])
                rightstr=midstr;
                midstr=[];
            end

            leftstr=[leftstr{:}];
            rightstr=[rightstr{:}];
            switch leftstr
                case 'int'
                    ocStruct.objective.integral.function.term=rightstr;
                    ocStruct.objective.integral.function.dependence=midstr;
                    ocStruct.objective.integral.function.term=rightstr;
                case 'expdisc'
                    ocStruct.objective.integral.discountrate=rightstr;
                    ocStruct.objective.integral.discountfactor.type=leftstr;
            end
        end

    case 'spacegeometry'
        for ii=1:numel(sectiontext)
            [leftstr,midstr,rightstr]=getpropstrings(sectiontext{ii});
            if isempty([rightstr{:}])
                rightstr=midstr;
                midstr=[];
            end

            leftstr=[leftstr{:}];
            %rightstr=[rightstr{:}];
            switch leftstr
                case 'intvl'
                    if length(rightstr)==2
                        ocStruct.spacegeometry.R1.interval.boundary.left.term=rightstr{1}(2:end);
                        ocStruct.spacegeometry.R1.interval.boundary.right.term=rightstr{2}(1:end-1);
                    end
            end
        end

    case 'boundarycondition'
        statename=retrieveppdemodelinformation(ocStruct,'statename');
        for ii=1:numel(sectiontext)
            [leftstr,midstr,rightstr]=getpropstrings(sectiontext{ii});
            if isempty([rightstr{:}])
                rightstr=midstr;
                midstr=[];
            end

            leftstr=[leftstr{:}];
            %rightstr=[rightstr{:}];
            switch leftstr
                case 'nm'
                    if ~isempty(midstr) && ~isempty(midstr{1})
                        switch midstr{1}
                            case 'intvl'
                                rightstr=[rightstr{:}];
                                RS=regexpi(rightstr,'(?<==)\S+','match');
                                RS=[RS{:}];
                                LS=regexpi(rightstr,'\S+(?==)','match');
                                LS=[LS{:}];
                                term=[LS '-(' RS ')'];
                                [normalderivativevariable,dum,idx1]=regexpi(term,'(?<=Dn.)\w+','match');
                                [normalderivativecoefficient,idx0,dum]=regexpi(term,'[+-]{1}\w+\*(?=Dn\.)|[+-]{1}\[\S+\]\*(?=Dn\.)','match');
                                if isempty(idx0)
                                    idx0=1;
                                    normalderivativecoefficient{1}='1';
                                end
                                term(idx0:idx1)=[];

                                fidx=find(strcmp(statename.value,normalderivativevariable));
                                ocStruct.constraint.derivative.state(fidx).boundarycondition.type='Neumann1D';
                                ocStruct.constraint.derivative.state(fidx).boundarycondition.left.coefficient=normalderivativecoefficient;
                                ocStruct.constraint.derivative.state(fidx).boundarycondition.left.term=term;
                                ocStruct.constraint.derivative.state(fidx).boundarycondition.right.coefficient=normalderivativecoefficient;
                                ocStruct.constraint.derivative.state(fidx).boundarycondition.right.term=term;

                        end
                    end
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

    case 'optimizationtype'
        ocStruct.optimizationtype=sectiontext{1};

    case 'controlconstraint'
        for ii=1:numel(sectiontext)
            [leftstr,midstr,rightstr]=getpropstrings(sectiontext{ii});
            midstr=[midstr{:}];
            leftstr=[leftstr{:}];
            rightstr=[rightstr{:}];
            ocStruct.constraint.function.control(ii).identifier=leftstr;
            ocStruct.constraint.function.control(ii).type=midstr;
            switch midstr
                case 'ineq'
                    fidx=strfind(rightstr,'>=');
                    if ~isempty(fidx)
                        leftidx=[1 fidx+2];
                        rightidx=[fidx-1 numel(rightstr)];
                        ocStruct.constraint.function.control(ii).term=[strtrim(rightstr(leftidx(1):rightidx(1))) '-(' strtrim(rightstr(leftidx(2):rightidx(2))) ')'];
                    end
                    fidx=strfind(rightstr,'<=');
                    if ~isempty(fidx)
                        leftidx=[1 fidx+2];
                        rightidx=[fidx-1 numel(rightstr)];
                        ocStruct.constraint.function.control(ii).term=[strtrim(rightstr(leftidx(2):rightidx(2))) '-(' strtrim(rightstr(leftidx(1):rightidx(1))) ')'];
                    end
            end
        end


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
