function v = ocmatfindsym(S,symkernel)
%FINDSYM Finds the symbolic variables in a symbolic expression or matrix.
%   FINDSYM(S), where S is a scalar or matrix sym, returns a string
%   containing all of the symbolic variables appearing in S. The
%   variables are returned in lexicographical order and are separated by
%   commas. If no symbolic variables are found, FINDSYM returns the
%   empty string.  The constants pi, i and j are not considered variables.
%
%   FINDSYM(S,N) returns the N symbolic variables closest to 'x'.
%
%   Examples:
%      findsym(alpha+a+b) returns
%       a, alpha, b
%
%      findsym(cos(alpha)*b*x1 + 14*y,2) returns
%       x1,y
%
%      findsym(y*(4+3*i) + 6*j) returns
%       y

%   Copyright 1993-2006 The MathWorks, Inc.
%   $Revision: 1.20.4.2 $  $Date: 2006/10/02 16:37:26 $

switch symkernel
    case 'maple'
        % Convert array to scalar.
        if iscell(S)
            sc = cell2vectorstring(S);
        else
            sc = char(S);
        end

        % Get the sorted list of all symbolic variables from Maple
        v = maple('indets', sc ,'symbol');

        % Return empty string if no symbolic variables found
        if isempty(v) || strcmp(v,'{}')
            v = '';
            return;
        end
        v=regexp(v(2:end-1),',[\ ]*','split');
        v(strcmp(v,'pi'))=[];
    case 'mupad'
        % Convert array to scalar.
        if verLessThan('symbolic','8')
            if iscell(S)
                sc = cell2vectorstring(S);
            else
                sc = char(S);
            end
            v=findsym(sym(sc));
            if isempty(v)
                v = '';
                return;
            end
            v=regexp(v,',[\ ]*','split');
            v(strcmp(v,'pi'))=[];
        else
            if iscell(S)
                sc = str2sym(cell2vectorstring(S));
            else
                sc = str2sym(S);
            end
            v=symvar(sc);
            if isempty(v)
                v = '';
                return;
            end
            v=cellfun(@char,sym2cell(v),'UniformOutput',false);
            v(strcmp(v,'pi'))=[];
        end
end