function varStruct=makestandardvariable4standardmodel(varStruct,transStruct,transformtype,varargin)
%
% transformtype: 0 ... symbolic, 1 ... array, 2 ... vectorized
% multline: if string is a matrix already decomposed in multiline format
% string is a cell array of strings
% if string contains ';' it is broken up into cells (standard matrix
% format)
% arcidentifier: is a string containing a digit identifying the specific
% arc. If arcidientifer is not empty the transformed variable names for
% this arc used for replacing. If arcidentifier is empty the standard
% transformation names are used for replacing

if nargin<=2
    transformtype=0;
end
fn=fieldnames(varStruct);
for ii=numel(fn):-1:1
    if ~strcmp(varStruct.(fn{ii}).type,'algterm')
        fn(ii)=[];
    end
end
if isempty(fn)
    return
end

switch transformtype
    case 0
        transformfieldname='name';
    case 1
        transformfieldname='arrayname';
    case 2
        transformfieldname='vectorizedname';
    otherwise
        ocmaterror('Transformtype unknown.')
end
for ii=1:numel(fn)
    %fprintf('Transform variabel ''%s'' (%d/%d)\n',fn{ii},ii,numel(fn))
    if isfield(varStruct.(fn{ii}),'string')
        for jj=1:numel(varStruct.(fn{ii}).string)
            if transformtype==2 && varStruct.(fn{ii}).vectorize
                % in case of vectorization values consisting only of parameter
                % values have to be adapted by using repmat
                if iscell(varStruct.(fn{ii}).string{jj})
                    numstr=numel(varStruct.(fn{ii}).string{jj});
                    novariable=ones(1,numstr);
                    nonumber=ones(1,numstr);
                    for kk=1:numel(transStruct.composedvar.username)
                        fres=regexp(varStruct.(fn{ii}).string{jj},['\<' transStruct.composedvar.username{kk} '\>']);
                        for ll=1:numstr
                            novariable(ll)=novariable(ll)&&isempty(fres{ll});
                        end
                        fres=regexp(varStruct.(fn{ii}).string{jj},'(^[0-9]+[;,\]\ \.]+$)|(^[\[\ ]+[0-9]+[;,\]\ \.]+$)|(^[A-z]+=[\[\ ]*[0-9]+[;,\]\ \.]+$)');
                        for ll=1:numstr
                            nonumber(ll)=nonumber(ll)&&~isempty(fres{ll});
                        end
                    end
                    fres=regexp(varStruct.(fn{ii}).string{jj},['\<' transStruct.independent.name '\>']);
                    for ll=1:numstr
                        novariable(ll)=novariable(ll)&&isempty(fres{ll});
                    end
                    if any(novariable)
                        fidx=find(novariable);
                        for idx=fidx
                            for mm=1:numel(transStruct.composedvar.parametername)
                                varStruct.(fn{ii}).string{jj}{idx}=regexprep(varStruct.(fn{ii}).string{jj}{idx},['\<' transStruct.composedvar.parametername{mm} '\>'],transStruct.composedvar.vectorizedparametername{mm});
                            end
                        end
                    end
                    if any(nonumber)
                        fidx=find(nonumber);
                        for idx=fidx
                            try
                                varStruct.(fn{ii}).string{jj}{idx}=regexprep(varStruct.(fn{ii}).string{jj}{idx},'([0-9]+)',['repmat($1,1,' transStruct.vectorizedsize.string ')']);
                            end
                        end
                    end
                else
                    novariable=1;
                    nonumber=1;
                    for kk=1:numel(transStruct.composedvar.username)
                        fres=regexp(varStruct.(fn{ii}).string{jj},['\<' transStruct.composedvar.username{kk} '\>']);
                        novariable=novariable&&isempty(fres);
                        fres=regexp(varStruct.(fn{ii}).string{jj},'(^[0-9]+[;,\]\ \.]+$)|(^[\[\ ]+[0-9]+[;,\]\ \.]+$)|(^[A-z]+=[\[\ ]*[0-9]+[;,\]\ \.]+$)');
                        nonumber=nonumber&&~isempty(fres);
                    end
                    fres=regexp(varStruct.(fn{ii}).string{jj},['\<' transStruct.independent.name '\>']);
                    novariable=novariable&&isempty(fres);
                    if novariable
                        for mm=1:numel(transStruct.composedvar.parametername)
                            varStruct.(fn{ii}).string{jj}=regexprep(varStruct.(fn{ii}).string{jj},['\<' transStruct.composedvar.parametername{mm} '\>'],transStruct.composedvar.vectorizedparametername{mm});
                        end
                    end
                    if nonumber
                        try
                            varStruct.(fn{ii}).string{jj}=regexprep(varStruct.(fn{ii}).string{jj},'([0-9]+)',['repmat($1,1,' transStruct.vectorizedsize.string ')']);
                        end
                    end
                end
            end
            if ~isempty(arcidentifier2field(varStruct.(fn{ii}).arcidentifier)) && varStruct.(fn{ii}).arcdependent>0
                if varStruct.(fn{ii}).arcdependent==1
                    transVar=transStruct.composedvar.(arcidentifier2field(varStruct.(fn{ii}).arcidentifier{jj}{1}));
                    numvar=numel(transVar.(transformfieldname));
                elseif varStruct.(fn{ii}).arcdependent==2
                    transVar=transStruct.composedvar.stdarc;
                    numvar=numel(transVar.(transformfieldname));
                elseif varStruct.(fn{ii}).arcdependent==3
                    transVar=transStruct.composedvar;
                    numvar=numel(transVar.(['standard' transformfieldname]));
                end
                for kk=1:numvar
                    if varStruct.(fn{ii}).arcdependent==3
                        varStruct.(fn{ii}).string{jj}=regexprep(varStruct.(fn{ii}).string{jj},['\<' strrep(strrep(strrep(transStruct.composedvar.username{kk},'(','\('),')','\)'),'_','\_') '\>'],transVar.(['standard' transformfieldname]){kk});
                    else
                        varStruct.(fn{ii}).string{jj}=regexprep(varStruct.(fn{ii}).string{jj},['\<' strrep(strrep(strrep(transStruct.composedvar.username{kk},'(','\('),')','\)'),'_','\_') '\>'],transVar.(transformfieldname){kk});
                    end
                end
            else
                if isfield(varStruct.(fn{ii}),'method') && strcmp(varStruct.(fn{ii}).method,'grad')
                    transVar=transStruct.locvar;
                    for kk=numel(transVar.(['standard' transformfieldname])):-1:1
                        varStruct.(fn{ii}).string{jj}=regexprep(varStruct.(fn{ii}).string{jj},['\<' strrep(strrep(strrep(transStruct.locvar.username{kk},'(','\('),')','\)'),'_','\_') '\>'],transVar.(['standard' transformfieldname]){kk});
                    end
                elseif ~isfield(varStruct,'method') || (isfield(varStruct,'method') && strcmp(varStruct.method,'bvp'))
                    transVar=transStruct.composedvar;
                    for kk=numel(transVar.(['standard' transformfieldname])):-1:1
                        varStruct.(fn{ii}).string{jj}=regexprep(varStruct.(fn{ii}).string{jj},['\<' strrep(strrep(strrep(transStruct.composedvar.username{kk},'(','\('),')','\)'),'_','\_') '\>'],transVar.(['standard' transformfieldname]){kk});
                    end
                end
%                 transVar=transStruct.composedvar;
%                 for kk=numel(transVar.(['standard' transformfieldname])):-1:1
%                     varStruct.(fn{ii}).string{jj}=regexprep(varStruct.(fn{ii}).string{jj},['\<' strrep(strrep(strrep(transStruct.composedvar.username{kk},'(','\('),')','\)'),'_','\_') '\>'],transVar.(['standard' transformfieldname]){kk});
%                 end
            end
        end
    end
end