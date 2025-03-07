function algebraic=algebraicodeterm2string(ocStruct,class,vectflag,varargin)
%
% ALGEBRAICTERM2STRING process the algebraic terms of the structure
% ocStruct describing the model into strings which are printed into the
% model files


algebraic.term='';
algebraic.type='';
algebraic.info=[];
arcfield='';
classtype='';

if isempty(ocStruct)
    return
end

if nargin==2
  vectflag=0;
end

if nargin>=4
  arcfield=varargin{1};
end

if nargin>=5
  classtype=varargin{2};
end
    
switch class

    case 'dynamics'
        dynamics=retrieveodemodelinformation(ocStruct,'dynamics','','');
        for ii=1:numel(dynamics.value)
            if ii==1 && numel(dynamics.value)>1
                algebraic.term{ii}=['[' dynamics.value{ii} '; ...'];
            elseif numel(dynamics.value)==1
                algebraic.term{ii}=dynamics.value{ii};
            elseif ii<numel(dynamics.value)
                algebraic.term{ii}=[dynamics.value{ii} '; ...'];
            else
                algebraic.term{ii}=[dynamics.value{ii} ']'];
            end
        end

    case 'exogenousfunctionterm'
        exogenousfunctionterm=retrieveodemodelinformation(ocStruct,'exogenousfunctionterm');
        exogenousfunctionnum=retrieveodemodelinformation(ocStruct,'exogenousfunctionnum');
        for ii=1:exogenousfunctionnum.value
            if ii==1
                if exogenousfunctionnum.value>1
                    algebraic.term{ii}=['out=[' exogenousfunctionterm.value{ii} '; ...'];
                else
                    algebraic.term{ii}=['out=' exogenousfunctionterm.value{ii} ';'];
                end
            elseif ii<exogenousfunctionnum.value
                algebraic.term{ii}=[exogenousfunctionterm.value{ii} '; ...'];
            else
                algebraic.term{ii}=[exogenousfunctionterm.value{ii} ']'];
            end
        end

    case 'symbolicexogenousfunctionname'
        exogenousfunctionname=retrieveodemodelinformation(ocStruct,'exogenousfunctionname');
        exogenousfunctionnum=retrieveodemodelinformation(ocStruct,'exogenousfunctionnum');
        symstr='mystr2sym';
        for ii=1:exogenousfunctionnum.value
            if ii==1 && exogenousfunctionnum.value>1
                algebraic.term{ii}=[symstr '([''[' exogenousfunctionname.value{ii} ';'' ...'];
            elseif exogenousfunctionnum.value==1
                algebraic.term{ii}=[symstr '(''' exogenousfunctionname.value{ii} ''')'];
            elseif ii<exogenousfunctionnum.value
                algebraic.term{ii}=['''' exogenousfunctionname.value{ii} ';'' ...'];
            else
                algebraic.term{ii}=['''' exogenousfunctionname.value{ii} ']''])'];
            end
        end

    case 'symbolicexogenousfunctionterm'
        exogenousfunctionterm=retrieveodemodelinformation(ocStruct,'exogenousfunctionterm','','');
        exogenousfunctionnum=retrieveodemodelinformation(ocStruct,'exogenousfunctionnum');
        if verLessThan('symbolic','8')
            symstr='mystr2sym';
        else
            symstr='mystr2sym';
        end
        for ii=1:exogenousfunctionnum.value
            if ii==1 && exogenousfunctionnum.value>1
                algebraic.term{ii}=[symstr '([''[' exogenousfunctionterm.value{ii} ';'' ...'];
            elseif exogenousfunctionnum.value==1
                algebraic.term{ii}=[symstr '(''' exogenousfunctionterm.value{ii} ''')'];
            elseif ii<exogenousfunctionnum.value
                algebraic.term{ii}=['''' exogenousfunctionterm.value{ii} ';'' ...'];
            else
                algebraic.term{ii}=['''' exogenousfunctionterm.value{ii} ']''])'];
            end
        end

    case 'symbolicdynamics'
        dynamics=retrieveodemodelinformation(ocStruct,'dynamics','','');
        if verLessThan('symbolic','8')
            symstr='mystr2sym';
        else
            symstr='mystr2sym';
        end
        for ii=1:numel(dynamics.value)
            if ii==1 && numel(dynamics.value)>1
                algebraic.term{ii}=[symstr '([''[' dynamics.value{ii} ';'' ...'];
            elseif numel(dynamics.value)==1
                algebraic.term{ii}=[symstr '(''' dynamics.value{ii} ''')'];
            elseif ii<numel(dynamics.value)
                algebraic.term{ii}=['''' dynamics.value{ii} ';'' ...'];
            else
                algebraic.term{ii}=['''' dynamics.value{ii} ']''])'];
            end
        end
        
    case 'dynamicsjacobian'
        dynamicsjacobian=retrieveodemodelinformation(ocStruct,'dynamicsjacobian',[],getsymkernel);
        numterm=numel(dynamicsjacobian.value);
        for ii=1:numterm
            if ii==1 && numel(dynamicsjacobian.value)>1
                algebraic.term{ii}=['[' dynamicsjacobian.value{ii} '; ...'];
            elseif numel(dynamicsjacobian.value)==1
                algebraic.term{ii}=['[' dynamicsjacobian.value{ii}];
            elseif ii<numterm
                algebraic.term{ii}=[dynamicsjacobian.value{ii} '; ...'];
            else
                algebraic.term{ii}=dynamicsjacobian.value{ii};
            end
        end
        algebraic.term{ii}=[algebraic.term{ii} ']'];
        
    case 'dynamicsparameterjacobian'
        dynamicsparameterjacobian=retrieveodemodelinformation(ocStruct,'dynamicsparameterjacobian',[],getsymkernel);
        numterm=numel(dynamicsparameterjacobian.value);
        for ii=1:numterm
            if ii==1 && numel(dynamicsparameterjacobian.value)>1
                algebraic.term{ii}=['[' dynamicsparameterjacobian.value{ii} '; ...'];
            elseif numel(dynamicsparameterjacobian.value)==1
                algebraic.term{ii}=['[' dynamicsparameterjacobian.value{ii}];
            elseif ii<numterm
                algebraic.term{ii}=[dynamicsparameterjacobian.value{ii} '; ...'];
            else
                algebraic.term{ii}=dynamicsparameterjacobian.value{ii};
            end
        end
        algebraic.term{ii}=[algebraic.term{ii} ']'];

    case 'symbolicdynamicsjacobian'
        symbolicdynamicsjacobian=retrieveodemodelinformation(ocStruct,'dynamicsjacobian',[],getsymkernel);
        numterm=numel(symbolicdynamicsjacobian.value);
        if verLessThan('symbolic','8')
            symstr='mystr2sym';
        else
            symstr='mystr2sym';
        end
        for ii=1:numterm
            if ii==1 && numel(symbolicdynamicsjacobian.value)>1
                algebraic.term{ii}=[symstr '([''[' symbolicdynamicsjacobian.value{ii} ';'' ...'];
            elseif numel(symbolicdynamicsjacobian.value)==1
                algebraic.term{ii}=[symstr '(''[' symbolicdynamicsjacobian.value{ii} ''''];
            elseif ii<numterm
                algebraic.term{ii}=['''' symbolicdynamicsjacobian.value{ii} ';'' ...'];
            else
                algebraic.term{ii}=['''' symbolicdynamicsjacobian.value{ii} ']'''];
            end
        end
        algebraic.term{ii}=[algebraic.term{ii} '])'];
        
    case 'symbolicdynamicsparameterjacobian'
        symbolicdynamicsparameterjacobian=retrieveodemodelinformation(ocStruct,'dynamicsparameterjacobian',[],getsymkernel);
        numterm=numel(symbolicdynamicsparameterjacobian.value);
        if verLessThan('symbolic','8')
            symstr='mystr2sym';
        else
            symstr='mystr2sym';
        end
        for ii=1:numterm
            if ii==1 && numel(symbolicdynamicsparameterjacobian.value)>1
                algebraic.term{ii}=[symstr '([''' symbolicdynamicsparameterjacobian.value{ii} ';'' ...'];
            elseif numel(symbolicdynamicsparameterjacobian.value)==1
                algebraic.term{ii}=[symstr '([''' symbolicdynamicsparameterjacobian.value{ii} ''''];
            elseif ii<numterm
                algebraic.term{ii}=['''' symbolicdynamicsparameterjacobian.value{ii} ';'' ...'];
            else
                algebraic.term{ii}=['''' symbolicdynamicsparameterjacobian.value{ii} ''''];
            end
        end
        algebraic.term{ii}=[algebraic.term{ii} '])'];

    case 'parametervariables'
        parametername=retrieveodemodelinformation(ocStruct,'parametername');
        algebraic.term=parametername.value;
end

if vectflag
    if iscell(algebraic.term)
        for ii=1:numel(algebraic.term)
            algebraic.term{ii}=simplevectorize(algebraic.term{ii});
        end
    elseif ischar(algebraic.term)
        algebraic.term=simplevectorize(algebraic.term);
    end
end

function v = simplevectorize(v)
%VECTORIZE Vectorize expression.
%   VECTORIZE(S), when S is a string expression, inserts a '.' before
%   any '^', '*' or '/' in S.  The result is a character string.
%
%   VECTORIZE(FUN), when FUN is an INLINE function object, vectorizes the
%   formula for FUN.  The result is the vectorized version of the INLINE
%   function.
%
%   See also INLINE/FORMULA, INLINE.

%   Copyright 1984-2002 The MathWorks, Inc. 
%   $Revision: 1.11 $  $Date: 2002/04/15 04:22:08 $

if isempty(v)
    v = [];
else
    for k = fliplr(find((v=='^') | (v=='*') | (v=='/')))
        v = [v(1:k-1) '.' v(k:end)];
    end
end