function filetext=processcommand(filetext,commandtype,varStruct,writefile)
%
% PROCESSCOMMAND text (cells of lines) including commands for the
% auto-generation of model files is processed and the processed text is
% returned. 

if nargin==2
    varStruct=[];
end
switch commandtype
    case 'IFCLAUSE'
        [sectionname,leftidx,rightidx]=findvariableintext(filetext,'IFCLAUSE');
        ifindex=[];
        ifcounter=0;
        endifcounter=0;
        for ii=1:numel(filetext)
            if strncmp(sectionname{ii},'!IF',3)
                ifcounter=ifcounter+1;
                ifindex(ifcounter,1)=ii;
            elseif strncmp(sectionname{ii},'!ENDIF',6)
                endifcounter=endifcounter+1;
                ifindex(endifcounter,2)=ii;
            end
        end
        ifindex=sortifclauses(ifindex);

        remline=ifindex(:).';
        % first process the IF clauses
        if ~isempty(ifindex)
            for ii=1:size(ifindex,1)
                actindex=ifindex(ii,1);
                if ~eval(sectionname{actindex}{1}((leftidx{actindex}+3):(rightidx{actindex}-1)))
                    remline=[remline ifindex(ii,1):ifindex(ii,2)];
                end
            end
        end
        filetext(remline)=[];
        fid_write=fopen(writefile,'w');
        if fid_write==-1
            ocmaterror(['Cannot open file ' writefile '.'])
        end
        linenum=numel(filetext);
        for ii=1:linenum
            fprintf(fid_write,'%s\n',filetext{ii});
        end
        fclose(fid_write);
        
    case 'CASESWITCH'
        sectionname=findvariableintext(filetext,'ADDARCCASE');
        varname=findvariableintext(filetext,'OCMATVARIABLE');

        fid_write=fopen(writefile,'w');
        if fid_write==-1
            ocmaterror(['Cannot open file ' writefile '.'])
        end
        linenum=numel(filetext);
        startnum=[];
        endnum=[];
        startcounter=0;
        endcounter=0;
        for ii=1:linenum
            if strcmp(sectionname{ii},'!STARTADDARCCASE!')
                startcounter=startcounter+1;
                startnum(startcounter)=ii;
            elseif strcmp(sectionname{ii},'!ENDADDARCCASE!')
                endcounter=endcounter+1;
                endnum(endcounter)=ii;
            end
        end
        if startcounter-endcounter
        end
        if ~isempty(startnum)
            % replace arcdependent variables
            counter=1;
            ii=0;
            while 1
                ii=ii+1;
                if startnum(counter)==ii
                    % search for ARCIDENTIFIER variable 
                    arcidentifiername=findvariableintext(filetext{ii+1},'OCMATVARIABLE');
                    arcidentifiername=char(arcidentifiername);
                    for arc=1:numel(varStruct.(arcidentifiername(2:end-1)).string)
                        for kk=(ii+1):(endnum(counter)-1)
                            newFileText=filetext{kk};
                            if ~isempty(varname{kk})
                                actvarname=varname{kk}{1}(2:end-1);
                                if ~varStruct.(actvarname).multline
                                    for jj=1:numel(varname{kk})
                                        actvarname=varname{kk}{jj}(2:end-1);
                                        newFileText=strrep(newFileText,[varname{kk}{jj}],varStruct.(actvarname).string{arc});
                                    end
                                    fprintf(fid_write,'%s\n',newFileText);
                                else
                                    for ll=1:numel(varStruct.(actvarname).string{arc})
                                        if ll==1
                                            tabulator=numel(regexp(newFileText,'(\t)'));
                                            newFileText=strrep(newFileText,[varname{kk}{1}],varStruct.(actvarname).string{arc}{ll});
                                        else
                                            fprintf(fid_write,repmat('\t',1,tabulator+1));
                                            newFileText=varStruct.(actvarname).string{arc}{ll};
                                        end
                                        fprintf(fid_write,'%s\n',newFileText);
                                    end
                                end
                            else
                                fprintf(fid_write,'%s\n',newFileText);
                            end
                            tabulator=0;
                        end
                    end
                    ii=endnum(counter);
                else
                    fprintf(fid_write,'%s\n',filetext{ii});
                end
                if ii==linenum
                    break
                end
            end
        else
            ii=0;
            while 1
                ii=ii+1;
                fprintf(fid_write,'%s\n',filetext{ii});
                if ii==linenum
                    break
                end
            end
        end
        fclose(fid_write);
        % read adapted file
        fid_read=fopen(writefile,'r');
        if fid_read==-1
            ocmaterror(['Cannot open file ' writefile '.'])
        end
        if verLessThan('matlab','7.6.0')
            filetext=textscan(fid_read,'%[^\n\r]','Whitespace',' \b','BufSize',1e8);
        else
            filetext=textscan(fid_read,'%[^\n\r]','Whitespace',' \b');
        end
        filetext=filetext{1};
        fclose(fid_read);


    case 'GENERALCASESWITCH'
        sectionname=findvariableintext(filetext,'ADDCASE');
        varname=findvariableintext(filetext,'OCMATVARIABLE');

        fid_write=fopen(writefile,'w');
        if fid_write==-1
            ocmaterror(['Cannot open file ' writefile '.'])
        end
        linenum=numel(filetext);
        startnum=[];
        endnum=[];
        startcounter=0;
        endcounter=0;
        for ii=1:linenum
            if strcmp(sectionname{ii},'!STARTADDCASE!')
                startcounter=startcounter+1;
                startnum(startcounter)=ii;
            elseif strcmp(sectionname{ii},'!ENDADDCASE!')
                endcounter=endcounter+1;
                endnum(endcounter)=ii;
            end
        end
        if startcounter-endcounter
        end
        if ~isempty(startnum)
            % replace arcdependent variables
            counter=1;
            ii=0;
            while 1
                ii=ii+1;
                if startnum(counter)==ii
                    % search for ARCIDENTIFIER variable 
                    arcidentifiername=findvariableintext(filetext{ii+1},'OCMATVARIABLE');
                    arcidentifiername=char(arcidentifiername);
                    for arc=1:numel(varStruct.(arcidentifiername(2:end-1)).string)
                        for kk=(ii+1):(endnum(counter)-1)
                            newFileText=filetext{kk};
                            if ~isempty(varname{kk})
                                actvarname=varname{kk}{1}(2:end-1);
                                if ~varStruct.(actvarname).multline
                                    for jj=1:numel(varname{kk})
                                        actvarname=varname{kk}{jj}(2:end-1);
                                        newFileText=strrep(newFileText,[varname{kk}{jj}],varStruct.(actvarname).string{arc});
                                    end
                                    fprintf(fid_write,'%s\n',newFileText);
                                else
                                    for ll=1:numel(varStruct.(actvarname).string{arc})
                                        if ll==1
                                            tabulator=numel(regexp(newFileText,'(\t)'));
                                            newFileText=strrep(newFileText,[varname{kk}{1}],varStruct.(actvarname).string{arc}{ll});
                                        else
                                            fprintf(fid_write,repmat('\t',1,tabulator+1));
                                            newFileText=varStruct.(actvarname).string{arc}{ll};
                                        end
                                        fprintf(fid_write,'%s\n',newFileText);
                                    end
                                end
                            else
                                fprintf(fid_write,'%s\n',newFileText);
                            end
                            tabulator=0;
                        end
                    end
                    ii=endnum(counter);
                else
                    fprintf(fid_write,'%s\n',filetext{ii});
                end
                if ii==linenum
                    break
                end
            end
        else
            ii=0;
            while 1
                ii=ii+1;
                fprintf(fid_write,'%s\n',filetext{ii});
                if ii==linenum
                    break
                end
            end
        end
        fclose(fid_write);
        % read adapted file
        fid_read=fopen(writefile,'r');
        if fid_read==-1
            ocmaterror(['Cannot open file ' writefile '.'])
        end
        if verLessThan('matlab','7.6.0')
            filetext=textscan(fid_read,'%[^\n\r]','Whitespace',' \b','BufSize',1e8);
        else
            filetext=textscan(fid_read,'%[^\n\r]','Whitespace',' \b');
        end
        filetext=filetext{1};
        fclose(fid_read);

    case 'FORNEXT'
        [sectionname,leftidx]=findvariableintext(filetext,'FORNEXT');
        if ~isempty([sectionname{:}])
            % process for next loop
            varname=findvariableintext(filetext,'OCMATVARIABLE');
            linenum=numel(filetext);
            forindex=zeros(sum([leftidx{:}]>0)/2,2);
            forcounter=0;
            endforcounter=0;
            for ii=1:linenum
                if strncmp(sectionname{ii},'!FOR',4)
                    forcounter=forcounter+1;
                    forindex(forcounter,1)=ii;
                elseif strncmp(sectionname{ii},'!ENDFOR',7)
                    endforcounter=endforcounter+1;
                    forindex(endforcounter,2)=ii;
                end
            end
            forindex(end+1,:)=inf;
            counter=1;
            ii=0;
            fid_write=fopen(writefile,'w');
            while 1
                ii=ii+1;
                if forindex(counter,1)==ii
                    forcountername=findvariableintext(filetext{ii},'FORNEXTCOUNTER');
                    forcountername=char(forcountername);
                    tmpout=regexp(filetext{ii},'=','split');
                    forcounteridx=str2num(tmpout{2}(1:end-1));
                    for jj=forcounteridx
                        for kk=(ii+1):(forindex(counter,2)-1)
                            newFileText=strrep(filetext{kk},forcountername,num2str(jj));
                            evalcommand=findvariableintext(newFileText,'EVALCOMMAND');
                            if ~isempty(evalcommand)
                                newFileText=strrep(newFileText,evalcommand{1},num2str(eval(regexprep(evalcommand{1}(2:end-1),'EVAL||\!||\!',''))));
                            end
                            if ~isempty(varname{kk})
                                actvarname=varname{kk}{1}(2:end-1);
                                if ~varStruct.(actvarname).multline
                                    for jj=1:numel(varname{kk})
                                        actvarname=varname{kk}{jj}(2:end-1);
                                        newFileText=strrep(newFileText,[varname{kk}{jj}],varStruct.(actvarname).string{arc});
                                    end
                                    fprintf(fid_write,'%s\n',newFileText);
                                else
                                    for ll=1:numel(varStruct.(actvarname).string{jj})
                                        if ll==1
                                            tabulator=numel(regexp(newFileText,'(\t)'));
                                            newFileText=strrep(newFileText,[varname{kk}{1}],varStruct.(actvarname).string{jj}{ll});
                                        else
                                            fprintf(fid_write,repmat('\t',1,tabulator+1));
                                            newFileText=varStruct.(actvarname).string{jj}{ll};
                                        end
                                        fprintf(fid_write,'%s\n',newFileText);
                                    end
                                end
                            else
                                fprintf(fid_write,'%s\n',newFileText);
                            end
                            %tabulator=2;
                        end
                    end
                    ii=forindex(counter,2);
                    counter=counter+1;
                else
                    fprintf(fid_write,'%s\n',filetext{ii});
                end
                if ii==linenum
                    fclose(fid_write);
                    break
                end
            end
        end
end

function ifindex=sortifclauses(ifindex)
% sort position of associated 'if' and 'end' clauses
if isempty(ifindex)
    return
end
if any(ifindex==0)
    ocmaterror('''If'' and ''end'' clause is not closing.')
end
ifclausepos=[ifindex(:,1)' ifindex(:,2)'];
ifclauseindex=[ones(1,size(ifindex,1)) zeros(1,size(ifindex,1))];
[ifclausepos sortidx]=sort(ifclausepos);
ifclauseindex=ifclauseindex(sortidx);

counter=size(ifindex,1);
while counter
    fidx=find(diff(ifclauseindex)==-1);
    for ii=length(fidx):-1:1
        ifindex(counter,:)=ifclausepos([fidx(ii) fidx(ii)+1]);
        counter=counter-1;
    end
    ifclausepos([fidx fidx+1])=[];
    ifclauseindex([fidx fidx+1])=[];
end