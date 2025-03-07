function S=ocmatsolve(eqns,vars,symkernel)
%OCMATSOLVE  Symbolic solution of algebraic equations.
%   OCMATSOLVE(eqns,varns)
%
%
%      S = solve('[x^2*y^2 - 2*x - 1,'x^2 - y^2 - 1]','[x,y]') returns
%      the solutions in a structure array.
%
%        S =
%          1x8 struct array with fields:
%          x
%          y
%

%   Copyright 1993-2007 The MathWorks, Inc.
%   $Revision: 1.21.4.3 $  $Date: 2007/06/14 05:26:44 $

% Set the Explicit solution (for equations of order 4 or less)
% option on in the Maple Workspace.

S=[];
if nargin==2
    symkernel='maple';
end

switch symkernel
    case 'maple'
        maple('_EnvExplicit := true;');
        
        % Collect input arguments together in either equation or variable lists.
        
        neqns = sum(commas(eqns)) + 1;
        nvars = sum(commas(vars)) + 1;
        
        if neqns ~= nvars
            ocmatmsg('%d equations in %d variables.',neqns,nvars)
        end
        eqns=strrep(strrep(eqns,'[','{'),']','}');
        vars=strrep(strrep(vars,'[','{'),']','}');
        % Seek analytic solution.
        
        % string format for 'solve' maple('solve','{eq1,eq2}','{x1,x2}')
        [R,stat] = maple('solve',eqns,vars);
        
        if ((isempty(R) || stat ~= 0) &&(nvars == neqns))
            % determine symbolic variables in equations
            v=regexp(removewhitespace(removecurledbrackets(maple('indets',eqns,'symbol'))),',','split');
            v(strcmp(v,'pi'))=[];
            total_vars=numel(v);
            if (nvars == total_vars)
                ocmatmsg('Unable to find analytic solution. Resort to numerical method ''maple/fsolve''\n')
                [R,stat] = maple('fsolve',eqns,vars);
            end
        end
        if stat
            ocmaterror('symbolic:solve:errmsg3',R)
        end
        
        % If still no solution, give up.
        if isempty(R)
            S='';
            return
        end
        
        % Expand any RootOf.
        
        while ~isempty(findstr(R,'RootOf'))
            k = min(findstr(R,'RootOf'));
            p = findstr(R,'{'); p = max(p(p<k));
            q = findstr(R,'}'); q = min(q(q>k));
            s = R(p:q);
            t = s(min(findstr(s,'RootOf'))+6:end);
            e = min(find(cumsum((t=='(')-(t==')'))==0));
            % RootOf with one argument, possibly an imbedded RootOf.
            checks = s;
            [s,stat] = maple('allvalues',s,'dependent');
            if isequal(checks,s)
                s = maple('evalf',s);
            end
            if isequal(checks,s)
                ocmatmsg('Unable to find closed form solution.')
                S='';
                return
            end
            R = [R(1:p-1) s R(q+1:end)];
        end
        
        % Parse the result.
        
        % Complete the structure.
        R=removewhitespace(R);
        % finde multiple solutions
        Rmult=regexp(R,'},{','split');
        
        % make structure skeleton
        v=regexp(removewhitespace(removecurledbrackets(vars)),',','split');
        for ii=1:nvars
            S(1).(v{ii})=[];
        end
        for ii=1:numel(Rmult)
            Rmult{ii}=removecurledbrackets(Rmult{ii});
            %Rmult{ii}=removeouterbrackets(Rmult{ii});
            solpart=regexp(Rmult{ii},',','split');
            for jj=1:numel(solpart)
                fidx=find(solpart{jj}=='=');
                if ~isempty(fidx)
                    S(ii).(solpart{jj}(1:fidx-1))=solpart{jj}(fidx+1:end);
                end
            end
        end
    case 'mupad'
        if verLessThan('symbolic','8')
            eqns=regexp(removeouterbrackets(eqns),',','split');
            idx=find(cellfun('isempty',regexp(eqns,'[*+-/^]')));
            for ii=idx
                eqns{ii}=[eqns{ii} '-0'];
            end
            vars=regexp(vars(2:end-1),',','split');
            nvars=length(vars);
            for ii=1:nvars
                S(1).(vars{ii})=[];
            end
            s=solve(eqns{:},vars{:});
            if isstruct(s)
                for ii=1:nvars
                    for jj=1:numel(s.(vars{ii}))
                        try
                            S(jj).(vars{ii})=char(simple(s.(vars{ii})(jj)));
                        catch
                            S(jj).(vars{ii})=char(simplify(s.(vars{ii})(jj)));
                        end
                    end
                end
            else
                for ii=1:numel(s)
                    try
                        S(ii).(vars{1})=char(simple(s(ii)));
                    catch
                        S(ii).(vars{1})=char(simplify(s(ii)));
                    end
                end
            end
        else
            eqns=regexp(removeouterbrackets(eqns),',','split');
            idx=find(cellfun('isempty',regexp(eqns,'[*+-/^]')));
            for ii=idx
                eqns{ii}=[eqns{ii} '-0'];
            end
            vars=regexp(vars(2:end-1),',','split');
            nvars=length(vars);
            for ii=1:nvars
                S(1).(vars{ii})=[];
            end
            varsName=vars;
            eqns=cellfun(@mystr2sym,eqns,'UniformOutput',false);
            vars=cellfun(@mystr2sym,vars,'UniformOutput',false);
            s=solve(eqns{:},vars{:},'IgnoreAnalyticConstraints',true,'ReturnConditions',false);
            if isstruct(s)
                for ii=1:nvars
                    for jj=1:numel(s.(varsName{ii}))
                        try
                            S(jj).(varsName{ii})=char(simple(s.(varsName{ii})(jj)));
                        catch
                            S(jj).(varsName{ii})=char(simplify(s.(varsName{ii})(jj)));
                        end
                    end
                end
            else
                for ii=1:numel(s)
                    try
                        S(ii).(varsName{1})=char(simple(s(ii)));
                    catch
                        S(ii).(varsName{1})=char(simplify(s(ii)));
                    end
                end
            end
        end
end
%-------------------------

function s = removewhitespace(s)
%TRIM  TRIM(s) deletes any white spaces.
s=regexprep(s,'\s','');

% function s = removeouterbrackets(s)
% %TRIM  TRIM(s) deletes any white spaces.
% s=regexprep(s,'([[)|(]])','');

function s = removecurledbrackets(s)
%TRIM  TRIM(s) deletes any white spaces.
s=regexprep(s,'{|}','');

%-------------------------

function c = commas(s)
%COMMAS  COMMAS(s) is true for commas not inside parentheses.
p = cumsum((s == '(') - (s == ')'));
c = (s == ',') & (p == 0);

