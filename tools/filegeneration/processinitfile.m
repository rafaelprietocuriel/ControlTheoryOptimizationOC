function [ocStruct,datafilename]=processinitfile(initfilename,opt)
%
% PROCESSINITFILE is the main file to read the initialization file
% describing the optimal control model.
%
% OCSTRUCT=PROCESSINITFILE(INITFILENAME) INITFILENAME is the name of
% initialization file without extension. A file 'INITFILENAME.EXT' has to
% exist in the initialization folder of 
% OCMat, where 'EXT' is one of the possible extensions 
%   'ocm'  .... (continuous) optimal control model
%   'docm'  .... discrete optimal control model
%   'odem'  .... dynamical system
% For the structure of the initialization file see examples in the
% 'initfiles' folder. OCSTRUCT is a structure containing the necessary
% information for a further processing and automatic file generation.
%
% OCSTRUCT=PROCESSINITFILE(INITFILENAME,OPT) OPT is an the ocoption
% structure.
% OPT.INIT.MessageDisplay='on'/'off': turning 'on'/'off' the more detailed
%       messaging about the processing of the initialization file 
% OPT.INIT.TestConsistency='on'/'off': turning 'on'/'off' the consistency
%       test for the initialization file.
%
% [OCSTRUCT,DATAFILENAME]=PROCESSINITFILE(...) DATAFILENAME is the full
% filename, where the structure OCSTRUCT is stored. The default path and
% filename are.
% ocmat/model/usermodel/[modelname]/data/[modelname]ModelDataStructure.mat

ocStruct=[];
if nargin==1
    opt=defaultocoptions;
end

[symkernel,symbolicinfo]=getsymkernel();

% get some options  
MessageDisplay=strcmp(getocoptions(opt,'INIT','MessageDisplay'),'on');
TestConsistency=strcmp(getocoptions(opt,'INIT','TestConsistency'),'on');

% 'initfilename' 
initfile=dir(fullocmatfile(getocmatfolder('initializationfile'),[initfilename '.*']));
if length(initfile)>1
    ocmaterror('Ambiguous file name')
elseif isempty(initfile)
    ocmaterror('No init file for ''%s'' does exist.',initfilename)
end
extension=char(regexp(initfile.name,'\.\w+\>','match'));   
if ~any(strcmp({'.ocm','.msm','.docm','.odem','.som','.dgm'},extension))
    ocmaterror('Extension ''%s'' of init file is unknown.',extension)
end
% remove dot of extension
extension(1)=[];

fread=fullocmatfile(getocmatfolder('initializationfile'),initfile.name);
if verLessThan('matlab','7.6.0')
    initfiletext=deblank(readfile2cell(fread,'%[^\n\r%]','CommentStyle','%','Whitespace',' \b\t','bufsize',4095*10));
else
    initfiletext=deblank(readfile2cell(fread,'%[^\n\r%]','CommentStyle','%','Whitespace',' \b\t'));
end
[modelname,modeltype,initsection,initfiletext,focsection,focfiletext]=separateintosections(initfiletext);
if isempty(modelname)
    modelname=initfilename;
else
    if ~(strcmp(modelname,initfilename))
        while 1
            answer=input(['Modelname ''' modelname ''' is different from initializationfile name ''' initfilename '''? Interrupt (y)/n: '],'s');
            if isempty(answer)
                % default value 'y'
                answer='y';
            end
            if strcmpi(answer,'n')
                break
            elseif strcmpi(answer,'y')
                return
            end
        end
    end
end
% depending on the type of model the specific files for reading the
% initialization file and deriving the necessary optimality conditions are
% called. The necessary optimality conditions can either be provided by the
% user in a section starting with the signal word FOC (see e.g.
% nonlincapitalaccumulation.ocm) or are derived using the symbolic toolbox
% of MATLAB (see e.g. mom.ocm, fireorwater.ocm)
initfunc=str2func([modeltype 'initialization']);
supplementfunc=str2func([modeltype 'supplement']);
readnocfunc=str2func([modeltype 'necessaryconditions']);
derivenocfunc=str2func([modeltype 'derivenecessaryconditions']);
testfunc=str2func([modeltype 'testconsistency']);

ocStruct.modeltype=modeltype;
ocStruct.modelname=modelname;
ocStruct.basicextension=extension(1:end-1);
ocStruct.initfilename=initfilename;

% test for syntax consistency of the initialization file
if TestConsistency
    b=testfunc(ocStruct,'maininitfile',symkernel,initsection,initfiletext);
    if ~b && MessageDisplay
        ocmatmsg('Initialization file passed the main syntax consistency test.\n\n')
    elseif b
        ocmatmsg('Initialization is aborted.\nInconsistency in the initialization file.\n')
        return
    end
    if ~isempty(focsection)
        b=testfunc(ocStruct,'focinitfile',symkernel,focsection,focfiletext);
        if ~b && MessageDisplay
            ocmatmsg('FOC part of the initialization file passed the first consistency test.\n\n')
        elseif b
            ocmatmsg('Initialization is aborted.\nInconsistency in FOC part of the initialization file.\n')
            return
        end
    end
end

% retrieve the model properties from the initialization file
ocStruct=initfunc(initsection,initfiletext,ocStruct);
if isempty(ocStruct)
    return
end
% complete missing entries for the main part in ocStruct by default values
ocStruct=supplementfunc(ocStruct,'maininitfile');

% test for consistency before the first order necessary conditions are
% derived, b=1 failed consistency test, b=0 passed consistency test
if TestConsistency
    b=testfunc(ocStruct,'basicmodelstructure',symkernel,initsection,initfiletext);
    if ~b && MessageDisplay
        ocmatmsg('Initialization file passed the consistency test for the derived model structure.\n\n')
    elseif b
        ocmatmsg('Initialization is aborted.\nInconsistency in the initialization file.\n\n')
        return
    end
end
% if focsection is not empty the first order necessary conditions (foc(s))
% are provided in the initialization file, otherwise the foc are derived
% using the symbolic toolbox (if present)
if ~isempty(focsection)
    % complete missing entries for the main part in ocStruct by default values
    ocStruct=supplementfunc(ocStruct,'focinitfile');    
    if MessageDisplay
        ocmatmsg('Start reading first order necessary optimality conditions, which are provided by the user.\n\n')
    end
    % the first order necessary optimality conditions are read out from the
    % initialization file
    ocStruct=readnocfunc(focsection,focfiletext,ocStruct);
else
    if isempty(symkernel) && MessageDisplay
        ocmatmsg('Symbolic toolbox is not activated.\nThe necessary conditions have to be provided in the initialization file.\n\n')
    else
        if MessageDisplay
            ocmatmsg('The necessary conditions are derived using the symbolic toolbox version %s.\n',symbolicinfo.Version)
        end
        % the first order necessary optimality conditions are derived using
        % the symbolic toolbox of MATLAB
        ocStruct=derivenocfunc(ocStruct,symkernel,opt);
        if  MessageDisplay
            ocmatmsg('The necessary conditions are successfully derived.\n\n')
        end
    end
end

% generate model folder
[modelfolder,datamodelfolder]=generatemodelfolder(modelname);
if MessageDisplay
    ocmatmsg('Model folder(s) are successfully created:\n%s\n%s\n\n',fullocmatfile(modelfolder),fullocmatfile(datamodelfolder))
end
% store the model data in the data folder of the model
datafilename=storemodeldata(ocStruct,datamodelfolder);
if MessageDisplay
    ocmatmsg('Data structure ''ocStruct'' is successfully stored in file:\n%s.\n',datafilename)
end
% add the model folder to the MATLAB path
searchpath=path;
modelfolder=fullocmatfile(getocmatfolder('specificmodel','',modelname));
if isempty(strfind(searchpath,[modelfolder pathsep]))
    while 1
        answer=input([strrep(modelfolder,filesep,'\\') ' is not on MATLAB path. Add it?  (y)/n: '],'s');
        if isempty(answer)
            % default value 'y'
            answer='y';
        end
        if strcmpi(answer,'n')
            break
        elseif strcmpi(answer,'y')
            addpath(modelfolder)
            savepath
            break
        end
    end
end
