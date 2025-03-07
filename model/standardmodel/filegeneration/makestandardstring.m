function newstring=makestandardstring(string,arcidentifier,transStruct,transformtype)
%
% transformtype: 0 ... array, 1 ... vectorized, 2 ... symbolic
% multline: if string is a matrix already decomposed in multiline format
% string is a cell array of strings
% if string contains ';' it is broken up into cells (standard matrix
% format) 
% arcidentifier: is a string containing a digit identifying the specific
% arc. If arcidientifer is not empty the transformed varaible names for
% this arc used for replacing. If arcidentifer is empty the standard
% transformation names are used for replacing

if nargin<=3
    transformtype=0;
end
string=makemultlinematrix(string);
if iscell(string)
    multline=1;
end
if transformtype==1
    if multline
        ocmatmsg('Matrix cannot be vectorized!')
    else
        string=vectorize(string);
    end
end
if multline
    if isempty(arcidentifier)
        for ll=1:numel(string)
            newstring.term{ll}=string{ll};
            for jj=1:numel(transStruct.composedvar.username)
                if transformtype==0
                    newstring.term{ll}=strrep(newstring.term{ll},transStruct.composedvar.username{jj},transStruct.composedvar.standardarrayname{jj});
                elseif transformtype==1
                    newstring.term{ll}=strrep(newstring.term{ll},transStruct.composedvar.username{jj},transStruct.composedvar.standardvectorizedname{jj});
                else
                    newstring.term{ll}=strrep(newstring.term{ll},transStruct.composedvar.username{jj},transStruct.composedvar.standardname{jj});
                end
            end
        end
    else
        for ii=1:numel(arcidentifier)
            arc=['arc' arcidentifier{ii}];
            for ll=1:numel(string)
                newstring.(arc).term{ll}=string{ll};
                for jj=1:numel(transStruct.composedvar.username)
                    if transformtype==0
                        newstring.(arc).term{ll}=strrep(newstring.(arc).term{ll},transStruct.composedvar.username{jj},transStruct.composedvar.(arc).arrayname{jj});
                    elseif transformtype==1
                        newstring.(arc).term{ll}=strrep(newstring.(arc).term{ll},transStruct.composedvar.username{jj},transStruct.composedvar.(arc).vectorizedname{jj});
                    else
                        newstring.(arc).term{ll}=strrep(newstring.(arc).term{ll},transStruct.composedvar.username{jj},transStruct.composedvar.(arc).name{jj});
                    end
                end
            end
        end
    end
else
    if isempty(arcidentifier)
        newstring.term=string;
        for jj=1:numel(transStruct.composedvar.username)
            if transformtype==0
                newstring.term=strrep(newstring.term,transStruct.composedvar.username{jj},transStruct.composedvar.standardarrayname{jj});
            elseif transformtype==1
                newstring.term=strrep(newstring.term,transStruct.composedvar.username{jj},transStruct.composedvar.standardvectorizedname{jj});
            else
                newstring.term=strrep(newstring.term,transStruct.composedvar.username{jj},transStruct.composedvar.standardname{jj});
            end
        end
    else
        for ii=1:numel(arcidentifier)
            arc=['arc' arcidentifier{ii}];
            newstring.(arc).term=string;
            for jj=1:numel(transStruct.composedvar.username)
                if transformtype==0
                    newstring.(arc).term=strrep(newstring.(arc).term,transStruct.composedvar.username{jj},transStruct.composedvar.(arc).arrayname{jj});
                elseif transformtype==1
                    newstring.(arc).term=strrep(newstring.(arc).term,transStruct.composedvar.username{jj},transStruct.composedvar.(arc).vectorizedname{jj});
                else
                    newstring.(arc).term=strrep(newstring.(arc).term,transStruct.composedvar.username{jj},transStruct.composedvar.(arc).name{jj});
                end
            end
        end
    end
end

function multlinestr=makemultlinematrix(str)

index=strfind(str,';');
if isempty(str)||isempty(index)||iscell(str)
    multlinestr=str;
    return
end
leftindex=[1 index+1];
rightindex=[index numel(str)];

for ii=1:numel(leftindex)-1
    multlinestr{ii}=[str(leftindex(ii):rightindex(ii)) ' ...'];
end
ii=ii+1;
multlinestr{ii}=str(leftindex(ii):rightindex(ii));