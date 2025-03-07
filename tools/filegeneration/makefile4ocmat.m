function modelfiles=makefile4ocmat(ocStruct,varargin)
%
% MAKEFILE4OCMAT generates the specific modelfiles from general template
% files.
%
% MODELFILES=MAKEFILE4OCMAT(OCSTRUCT) OCSTRUCT is the model structure
% generated and stored during the processing of the initialization file.
% The specific model files, e.g. the file fore the canonical system,
% Hamiltonian, etc., are build from general template files. These template
% files can be found in 'ocmat\model\standardmodel\templatefiles'. On basis
% of these template files the user can create her/his own template files.
%
% MODELFILES=MAKEFILE4OCMAT(OCSTRUCT,FN) FN is a character variable or cell
% of character variables providing the filenames that are generated.
%
% MODELFILES=MAKEFILE4OCMAT(OCSTRUCT,FN,OPT) specific options OPT can be
% provided.
%
% MODELFILES=MAKEFILE4OCMAT(OCSTRUCT,FN,OPT,'EXCLUDE',EXCFN) the files
% EXCFN can be excluded from the file generation process.
%
% In the m-file GETMODELFILENAMES the filenames and some properties are
% defined that are then actually created. MODELFILES is a structure array
% that is returned by GETMODELFILENAMES.
% MODELFILES.NAME: basic name of a specific file (e.g. 'OptimalControl')
%           .TRANSFORMTYPE: determines if the terms are considered
%               vectorized/symbolic form
%           .FILENAME: full name of the specific file (e.g.
%               'F:\Matlab\ocmat\model\usermodel\out\momOptimalControl.m')

modelfiles=[];
filename='';
opt=[];
onlyreturnnames=0;
if isempty(ocStruct)
    return
end
if nargin>=2
    filename=varargin{1};
end
if nargin>=3
    opt=varargin{2};
end
if nargin>=4
    onlyreturnnames=varargin{3};
end
if isempty(opt)
    opt=defaultocoptions;
end
if isempty(onlyreturnnames)
    onlyreturnnames=0;
end
JacobianCalc=getocoptions(opt,'INIT','Jacobian');
HessianCalc=getocoptions(opt,'INIT','Hessian');
if strcmp(modeltype(ocStruct),'standardmodel')
    if implicitcontrols(ocStruct) && (strcmpi(JacobianCalc,'explicit') || strcmpi(HessianCalc,'explicit'))
        fprintf('\nThe symbolic derivation of the Jacobian/Hessian in case of implicit controls can take a considerably long time.\n')
        fprintf('To calculate the Jacobian/Hessian numerically set\n\n')
        fprintf('opt=setocoptions(''INIT'',''Jacobian'',''numerical'',''Hessian'',''numerical'');\n\n')
        fprintf('and restart the file generating process:\n\n')
        fprintf('modelfiles=makefile4ocmat(ocStruct,[],opt);\n\n')
        answer=input('Interrupt the processing (y)/n: ','s');
        while 1
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

%     if strcmpi(JacobianCalc,'numerical')
%         excludeidx=find(strcmpi(varargin,'exclude'));
%         includeidx=find(strcmpi(varargin,'include'));
%         if ~isempty(includeidx)
%             includefiles=varargin{includeidx+1};
%             if ischar(includefiles)
%                 includefiles={includefiles};
%             end
%             varargin{includeidx+1}=[includefiles {'CanonicalSystemNumericJacobian','EquilibriumEquationNumericJacobian','EquilibriumEquationNumericParameterJacobian'}];
%         else
%             varargin{end+1}='include';
%             varargin{end+1}={'CanonicalSystemNumericalParameterJacobian','CanonicalSystemNumericalJacobian','EquilibriumEquationNumericalJacobian','EquilibriumEquationNumericalParameterJacobian'};
%         end
%         if ~isempty(excludeidx)
%             excludefiles=varargin{excludeidx+1};
%             if ischar(excludefiles)
%                 excludefiles={excludefiles};
%             end
%             varargin{excludeidx+1}=[excludefiles {'SymbolicCanonicalSystemJacobian','EquilibriumEquationJacobian','EquilibriumEquationParameterJacobian', ...
%                 'CanonicalSystemImplicitJacobian','CanonicalSystemJacobian','CanonicalSystemParameterJacobian','CanonicalSystemImplicitNumericalJacobian', ...
%                 'CanonicalSystemImplicitNumericalParameterJacobian'}];
%         else
%             varargin{end+1}='exclude';
%             varargin{end+1}={'SymbolicCanonicalSystemJacobian','EquilibriumEquationJacobian','EquilibriumEquationParameterJacobian', ...
%                 'CanonicalSystemImplicitJacobian','CanonicalSystemJacobian','CanonicalSystemParameterJacobian', ...
%                 'CanonicalSystemImplicitNumericalJacobian','CanonicalSystemImplicitNumericalParameterJacobian'};
%         end
%     end
    if strcmpi(JacobianCalc,'implicit')
        excludeidx=find(strcmpi(varargin,'exclude'));
        includeidx=find(strcmpi(varargin,'include'));
        if ~isempty(includeidx)
            includefiles=varargin{includeidx+1};
            if ischar(includefiles)
                includefiles={includefiles};
            end
            varargin{includeidx+1}=[includefiles {'ImplicitGXNumerical','ImplicitGUNumerical','ImplicitGPNumerical'}];
        else
            varargin{end+1}='include';
            varargin{end+1}={'ImplicitGXNumerical','ImplicitGUNumerical','ImplicitGPNumerical'};
        end
        if ~isempty(excludeidx)
            excludefiles=varargin{excludeidx+1};
            if ischar(excludefiles)
                excludefiles={excludefiles};
            end
            varargin{excludeidx+1}=[excludefiles {'ImplicitGX','ImplicitGU','ImplicitGP'}];
        else
            varargin{end+1}='exclude';
            varargin{end+1}={'ImplicitGX','ImplicitGU','ImplicitGP'};
        end
    end

    if strcmpi(HessianCalc,'numerical')
        excludeidx=find(strcmpi(varargin,'exclude'));
        if ~isempty(excludeidx)
            excludefiles=varargin{excludeidx+1};
            if ischar(excludefiles)
                excludefiles={excludefiles};
            end
            varargin{excludeidx+1}=[excludefiles {'CanonicalSystemHessian','CanonicalSystemTotalHessian'}];
        else
            varargin{end+1}='exclude';
            varargin{end+1}={'CanonicalSystemHessian','CanonicalSystemTotalHessian'};
        end
    end
end
modelfiles=getmodelfilenames(ocStruct,varargin{:});
if ~isempty(filename)
    if ~iscell(filename)
        if strcmp(filename(end),'*')
            modelfiles=modelfiles(strncmp({modelfiles.name},filename,length(filename)-1));
        else
            modelfiles=modelfiles(strcmp({modelfiles.name},filename));
        end
    else
        keepidx=zeros(1,length(filename));
        for ii=1:length(filename)

            if strcmp(filename(end),'*')
                idx=find(strncmp({modelfiles.name},filename{ii},length(filename{ii})-1));
            else
                idx=find(strcmp({modelfiles.name},filename{ii}));
            end
            if ~isempty(idx)
                keepidx(ii)=idx;
            end
        end
        keepidx(keepidx==0)=[];
        modelfiles=modelfiles(keepidx);
    end
    if isempty(modelfiles)
        return
    end
end
if onlyreturnnames
    return
end
numfile=numel(modelfiles);

% the transformation structure consists of translation rules from user
% specifc notation to standard notation in the auto-generated files
transStruct=[];
if isfield(ocStruct,'optimization')
    for ii=1:size(ocStruct.optimization.method,1)
        switch deblank(ocStruct.optimization.method(ii,:))
            case 'bvp'
                varname={'independent','parameter','composedvar','vectorizedsize','time','space'};
                for jj=1:numel(varname)
                    transStruct=generatetransformationstruct(ocStruct,transStruct,varname{jj});
                end
            case 'dae'
                varname={'independent','parameter','composedvar','vectorizedsize','time','space'};
                for jj=1:numel(varname)
                    transStruct=generatetransformationstruct(ocStruct,transStruct,varname{jj});
                end
            case 'grad'
                varname={'independent','parameter','composedvar','vectorizedsize','time','space'};
                for jj=1:numel(varname)
                    transStruct=generatetransformationstruct(ocStruct,transStruct,varname{jj});
                end
                varname={'locvar','vectorizedsize','age','space'};
                for jj=1:numel(varname)
                    transStruct=generatetransformationstruct4grad(ocStruct,transStruct,varname{jj});
                end
            case {'zero','statgrad','semismooth'}
                varname={'variable','composedvar'};
                for jj=1:numel(varname)
                    transStruct=generatetransformationstruct4staticoptmodel(ocStruct,transStruct,varname{jj});
                end
        end
    end
else
    varname={'independent','parameter','composedvar','vectorizedsize','time','space'};
    for jj=1:numel(varname)
        transStruct=generatetransformationstruct(ocStruct,transStruct,varname{jj});
    end
end
% in a first pass the template file is searched for variables (these are of
% the form $VARIABLE$) and a structure consisting of theses variables and
% replacements is generated

% generate option dependent variables
varStruct=generatevariablestruct4options(opt);
%varStruct=[];
% generate model dependent variables
for ii=1:numfile
    varStruct=findvariable(modelfiles(ii).name,ocStruct,varStruct);
end

% the variable terms retrieved in the previous pass and stored in varStruct
% are transformed into standard notation, 0 replaces user notation to
% standard notation in symbolic form, 1 replaces to coordinate form and 2
% replaces to vectorized form.
varTransStruct0=makestandardvariable(ocStruct,varStruct,transStruct,0);
varTransStruct1=makestandardvariable(ocStruct,varStruct,transStruct,1);
varTransStruct2=makestandardvariable(ocStruct,varStruct,transStruct,2);

% the m-files are generated from the template files using the model
% specific values.
for ii=1:numfile
    switch modelfiles(ii).transformtype
        case 0
            modelfiles(ii).filename=writetemplate2mfile(modelfiles(ii).name,ocStruct,varTransStruct0);
        case 1
            modelfiles(ii).filename=writetemplate2mfile(modelfiles(ii).name,ocStruct,varTransStruct1);
        case 2
            modelfiles(ii).filename=writetemplate2mfile(modelfiles(ii).name,ocStruct,varTransStruct2);
        case 3
            modelfiles(ii).filename=writetemplate2mfile(modelfiles(ii).name,ocStruct,varStruct);
    end
end
