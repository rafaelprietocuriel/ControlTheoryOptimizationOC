function structurebranching(ocStruct,valueflag)
% function fn_structdisp Xname
% function fn_structdisp(X)
%---
% Recursively display the content of a structure and its sub-structures
%
% Input:
% - Xname/X     one can give as argument either the structure to display or
%               or a string (the name in the current workspace of the
%               structure to display)
%
% A few parameters can be adjusted inside the m file to determine when
% arrays and cell should be displayed completely or not

% License text
% Copyright (c) 2005, Thomas Deneux
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without 
% modification, are permitted provided that the following conditions are 
% met:
% 
%     * Redistributions of source code must retain the above copyright 
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright 
%       notice, this list of conditions and the following disclaimer in 
%       the documentation and/or other materials provided with the distribution
%       
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
% POSSIBILITY OF SUCH DAMAGE.

% Thomas Deneux
% Copyright 2005-2012

if nargin==1
    valueflag=1;
end
structName = inputname(1);

outputfile=fullocmatfile(getocmatfolder('userdata',ocStruct.modeltype,ocStruct.modelname),['recursiveocstruct' num2str(valueflag) '.txt']);
try
    fclose(fopen(outputfile,'w+'));
end
diary(outputfile);
rec_structdisp(structName,ocStruct,valueflag)
diary off

%---------------------------------
function rec_structdisp(Xname,X,valueflag)
%---

%-- PARAMETERS (Edit this) --%

ARRAYMAXROWS = 10;
ARRAYMAXCOLS = 50;
ARRAYMAXELEMS = 500;
CELLMAXROWS = 10;
CELLMAXCOLS = 50;
CELLMAXELEMS = 500;
%CELLRECURSIVE = false;
CELLRECURSIVE = true;

%----- PARAMETERS END -------%
disp([Xname ':'])
if valueflag
    disp(X)
end
%fprintf('\b')

if isstruct(X) || isobject(X)
    F = fieldnames(X);
    nsub = length(F);
    Y = cell(1,nsub);
    subnames = cell(1,nsub);
    for i=1:nsub
        f = F{i};
        Y{i} = X.(f);
        subnames{i} = [Xname '.' f];
    end
elseif CELLRECURSIVE && iscell(X)
    nsub = numel(X);
    s = size(X);
    Y = X(:);
    subnames = cell(1,nsub);
    for i=1:nsub
        inds = s;
        globind = i-1;
        for k=1:length(s)
            inds(k) = 1+mod(globind,s(k));
            globind = floor(globind/s(k));
        end
        subnames{i} = [Xname '{' num2str(inds,'%i,')];
        subnames{i}(end) = '}';
    end
else
    return
end

for i=1:nsub
    a = Y{i};
    if isstruct(a) || isobject(a)
        if length(a)==1
            rec_structdisp(subnames{i},a,valueflag)
        else
            for k=1:length(a)
                rec_structdisp([subnames{i} '(' num2str(k) ')'],a(k),valueflag)
            end
        end
    elseif iscell(a)
        if size(a,1)<=CELLMAXROWS && size(a,2)<=CELLMAXCOLS && numel(a)<=CELLMAXELEMS
            rec_structdisp(subnames{i},a,valueflag)
        end
    elseif ischar(a)
        disp([subnames{i} ':'])
        if valueflag
            disp(a)
        end
    elseif size(a,1)<=ARRAYMAXROWS && size(a,2)<=ARRAYMAXCOLS && numel(a)<=ARRAYMAXELEMS
        disp([subnames{i} ':'])
        if valueflag
            disp(a)
        end
    end
end
