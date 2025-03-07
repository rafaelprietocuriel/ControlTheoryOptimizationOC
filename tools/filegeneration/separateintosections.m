function [modelname,modeltype,initsection,inittext,focsection,foctext]=separateintosections(inittext)

initsection=[];
focsection=[];
modelname='';
modeltype='';

if isempty(inittext)
    return
end
if ~iscell(inittext)
    ocmatmsg('Input argument is not a cell.')
    return
end
numline=numel(inittext);
% determine type of problem
fidx=find(strcmpi(inittext,'type'));
if ~isempty(fidx) && numel(fidx)==1 && numline>fidx
    modeltype=inittext{fidx+1};
end
if isempty(modeltype)
    ocmatmsg('Section ''Type'' is missing.')
    return
end
% determine model name
fidx=find(strcmpi(inittext,'modelname'));
if ~isempty(fidx) && numel(fidx)==1 && numline>fidx
    modelname=inittext{fidx+1};
end

% determine if first order necessary optimality conditions are provided
fidx=find(strcmpi(inittext,'foc'));
foctext=inittext(fidx+1:numline);
inittext(fidx:numline)=[];

initsection=text2structure(inittext,feval(str2func([modeltype 'properties']),'init'));
if isempty(initsection)
    return
end
if numel(foctext)
    focsection=text2structure(foctext,feval(str2func([modeltype 'properties']),'foc'));
    if isempty(focsection)
        initsection=[];
        return
    end
end

function sectionstruct=text2structure(celltext,propnames)

sectionstruct=[];
numprop=numel(propnames);
sortsec=zeros(numprop,1);
proporder=zeros(numprop,1);
sectionfieldname=cell(numprop,1);
for ii=1:numprop
    fidx=find(strcmpi(celltext,propnames{ii}));
    if ~isempty(fidx)
        if numel(fidx)>1
            ocmatmsg('The section ''%s'' is multiple defined\nInitialization process aborted',propnames{ii})
            return
        end
        sortsec(ii)=fidx;
        proporder(ii)=ii;
        sectionfieldname{ii}=propnames{ii};
    end
end
removeidx=~sortsec;
sortsec(removeidx)=[];
proporder(removeidx)=[];
sectionfieldname(removeidx)=[];
[sortsec sortidx]=sort(sortsec);
sortsec=[sortsec+1 [sortsec(2:end)-1;numel(celltext)]];
sectionfieldname=sectionfieldname(sortidx);
for ii=1:size(sortsec,1)
    sectionstruct.(sectionfieldname{ii})=sortsec(ii,:);
end
sectionstruct=orderfields(sectionstruct,propnames(proporder));
for ii=1:size(sortsec,1)
    if diff(sectionstruct.(sectionfieldname{ii}))<0
        sectionstruct=rmfield(sectionstruct,sectionfieldname{ii});
    end
end
