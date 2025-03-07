function ocStruct=odemodelinitialization(initsection,initfiletext,ocStruct)
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
                otherwise
                    ocStruct.variable.(leftstr).name=rightstr;
                    ocStruct.variable.(leftstr).num=numel(rightstr);
            end
        end
        
    case 'dynamics'
        odecounter=0;
        statevarname=ocStruct.variable.state.name;
        for ii=1:numel(sectiontext)
            [leftstr rightstr]=getpropstrings(sectiontext{ii});
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
