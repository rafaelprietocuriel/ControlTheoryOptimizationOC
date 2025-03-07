function ocStruct=standardmodelnecessaryconditions(focsection,foctext,ocStruct)
%
% reading out the necessary optimality conditions from the initialization
% file

if ~isfield(ocStruct,'foc')
    ocStruct.foc=[];
end
propname=fieldnames(focsection);
for ii=1:numel(propname)
    try
        ocStruct=propertyhandling(ocStruct,foctext(focsection.(propname{ii})(1):focsection.(propname{ii})(2)),propname{ii});
    catch
        ocStruct=[];
        ocmatmsg('Initialization process aborted.\nProblems in reading FOC section ''%s''.\n',propname{ii})
        return
    end
end


function ocStruct=propertyhandling(ocStruct,sectiontext,propname)

if isfield(ocStruct.foc,'abbreviation')
    abbreviation=ocStruct.foc.abbreviation;
else
    abbreviation=[];
end
switch propname
    case 'abbreviation'
        for ii=1:numel(sectiontext)
            [leftstr middlestr rightstr]=getpropstrings(sectiontext{ii});
            leftstr=[leftstr{:}];
            middlestr=[middlestr{:}];
            rightstr=[rightstr{:}];
            if strcmp(leftstr,'abb')
                if isempty(rightstr)
                    % if no abbreviation term is provided it is assumed
                    % that it is the identifier of a constraint
                    ccidentifier=retrievemodelinformation(ocStruct,'inequalitycontrolconstraintidentifier');
                    ccterm=retrievemodelinformation(ocStruct,'inequalitycontrolconstraint');
                    abbreviation.name{ii}=middlestr;
                    fidx=find(strcmp(middlestr,ccidentifier.value));
                    abbreviation.replacestr{ii}=ccterm.value{fidx};
                else
                    abbreviation.name{ii}=middlestr;
                    abbreviation.replacestr{ii}=rightstr;
                end
            end
        end
        ocStruct.foc.abbreviation=abbreviation;
    case 'pontryaginfunction'

        arcidentifier=ocStruct.arc.identifier;
        for ii=1:numel(sectiontext)
            [leftstr rightstr]=getpropstrings(sectiontext{ii});
            leftstr=[leftstr{:}];
            rightstr=[rightstr{:}];
            ocStruct.pontryaginfunction.identifier=leftstr;
            rightstr=removebrackets(replaceabbreviation(rightstr,abbreviation));
            ocStruct.pontryaginfunction.term=rightstr;
        end
        ocStruct.pontryaginfunction.arcidentifier=arcidentifier;
        
    case {'pontryaginfunctionderivativeorder1','pontryaginfunctionderivativeorder2'}
        linecounter=1;
        while linecounter<numel(sectiontext)
            [leftstr rightstr]=getpropstrings(sectiontext{linecounter});
            leftstr=[leftstr{:}];
            [multlinestr numlines]=findmultlinematrix([rightstr{:}],sectiontext,linecounter);
            linecounter=linecounter+numlines;
            rightstr=matrix2cell(removebrackets(replaceabbreviation(multlinestr,abbreviation)));
            ocStruct.pontryaginfunction.derivative.(leftstr).identifier=['D' ocStruct.pontryaginfunction.identifier leftstr];
            ocStruct.pontryaginfunction.derivative.(leftstr).term=rightstr;
        end

    case 'maximizingcondition'
        ccnum=retrievemodelinformation(ocStruct,'inequalitycontrolconstraintnum');
        for ii=1:numel(sectiontext)
            [leftstr middlestr rightstr]=getpropstrings(sectiontext{ii});
            for jj=1:numel(leftstr) % arc
                arc=arcidentifier2field(leftstr{jj});
                ocStruct.variable.control.(arc).arcidentifier=leftstr(jj);
                if ~isempty(ccnum.value)
                        ocStruct.foc.value.lagrangemultcc.(arc).arcidentifier=leftstr(jj);
                end
                for jj=1:numel(middlestr) % control
                    fidx=find(strcmp(middlestr{jj},ocStruct.variable.control.name));
                    if ~isempty(fidx)
                        rightstr=replaceabbreviation(rightstr{:},abbreviation);
                        ocStruct.foc.value.control.(arc).term{fidx}=rightstr;
                    end
                    fidx=find(strcmp(middlestr{jj},ocStruct.variable.lagrangemultcc.name));
                    if ~isempty(fidx)
                        rightstr=replaceabbreviation(rightstr{:},abbreviation);
                        ocStruct.foc.value.lagrangemultcc.(arc).term{fidx}=rightstr;
                    end
                end
            end
        end
    case 'adjointsystem'
        linecounter=1;
        while linecounter<numel(sectiontext)
            [leftstr middlestr rightstr]=getpropstrings(sectiontext{linecounter});
            [multlinestr numlines]=findmultlinematrix([rightstr{:}],sectiontext,linecounter);
            linecounter=linecounter+numlines;
            rightstr=vector2cell(removebrackets(replaceabbreviation(multlinestr,abbreviation)));
            if strcmp([leftstr{:}],'*')
                leftstr=ocStruct.arc.identifier;
            elseif ~isempty(strfind([leftstr{:}],':'))
                arcid=[eval(leftstr{:})];
                for ii=1:numel(arcid)
                    leftstr{ii}=num2str(arcid(ii));
                end
            end
            arc=arcidentifier2field(leftstr);
            middlestr=[middlestr{:}];
            ocStruct.foc.adjointsystem.equation.(arc).(middlestr).term=rightstr;
            ocStruct.foc.adjointsystem.equation.(arc).(middlestr).num=numel(rightstr);
        end

    case {'canonicalsystemderivativeorder1','canonicalsystemderivativeorder2'}
        linecounter=1;
        while linecounter<numel(sectiontext)
            [leftstr middlestr rightstr]=getpropstrings(sectiontext{linecounter});
            if strcmp([leftstr{:}],'*')
                leftstr=ocStruct.arc.identifier;
            end
            if strcmp([leftstr{:}],'*')
                leftstr=ocStruct.arc.identifier;
            end
            arc=arcidentifier2field(leftstr);
            middlestr=[middlestr{:}];
            [multlinestr numlines]=findmultlinematrix([rightstr{:}],sectiontext,linecounter);
            linecounter=linecounter+numlines;
            rightstr=matrix2cell(removebrackets(replaceabbreviation(multlinestr,abbreviation)));
            ocStruct.foc.canonicalsystem.derivative.(arc).(middlestr).term=rightstr;
        end


end

function varargout=getpropstrings(inputstr)
varargout=cell(nargout,1);
tmpcell=regexp(inputstr,'::','split');
if ~isempty(tmpcell)
    varargout(1:numel(tmpcell))=deblank(tmpcell);
else
    varargout{1}='';
    varargout{2}=inputstr;
end

% if output argument is comma seperated transform it into cell array beside
% the outmost right argument, which possibly is a mathematical term in
% MATLAB notation
for ii=1:numel(varargout)-1
    if ~isempty(strfind(varargout{ii},','))
        varargout{ii}=regexp(varargout{ii},',','split');
    else
        varargout{ii}={varargout{ii}};
    end
end
varargout{end}={varargout{end}};

function str=replaceabbreviation(str,abbreviation)

if isempty(abbreviation)
    return
end

for ii=1:numel(abbreviation.name)
    if size(abbreviation.replacestr{ii},1)>1
        % if replacestr is a multiline matrix the variable str is replaced
        % by this matrix
        str=abbreviation.replacestr{ii};
        return
    else
        str=strrep(str,abbreviation.name{ii},abbreviation.replacestr{ii});
    end
end

function [multlinestr,counter]=findmultlinematrix(str,text,pos)

counter=1;
if isempty(str)||isempty(regexp(str,'(;\ \.\.\.)$'))
    multlinestr=str;
    return
end

multlinestr{1}=str;
while 1
    pos=pos+1;
    initline=text{pos};
    counter=counter+1;
    multlinestr{counter}=initline;
    if strcmp(initline(end),']')
        return
    end
end

function str=removebrackets(str)

if isempty(str)
    return
end
str=regexprep(str,'\[|\]|(;\ \.\.\.)','');

function cellmatrix=matrix2cell(str)

if iscell(str) && numel(str)>1
    cellmatrix=str;
    return
end
cellmatrix=regexp(str,';','split');

function cellvector=vector2cell(str)

if iscell(str) && numel(str)>1
    cellvector=str;
    return
end
cellvector=regexp(str,'[;,\ ]','split');
