function algebraic=algebraicstandarddiffmodelterm2string(ocStruct,class,vectflag,varargin)
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
    case 'controlvalue'
        optimalvalue=retrievediffmodelinformation(ocStruct,'optimalvalue',field2arcidentifier(arcfield));
        controlnum=retrievediffmodelinformation(ocStruct,'controlnum');
        controlnametp1=retrievediffmodelinformation(ocStruct,'controlnametp1');
        for ii=1:controlnum.value
            if ii==1
                if controlnum.value>1
                    algebraic.term{ii}=['[' optimalvalue.value.(controlnametp1.value{ii}) '; ...'];
                else
                    algebraic.term{ii}=optimalvalue.value.(controlnametp1.value{ii});
                end
            elseif ii<controlnum.value
                algebraic.term{ii}=[optimalvalue.value.(controlnametp1.value{ii}) '; ...'];
            else
                algebraic.term{ii}=[optimalvalue.value.(controlnametp1.value{ii})];
            end
        end
        if controlnum.value>1
            algebraic.term{ii}=[algebraic.term{ii} ']'];
        end

    case 'lagrangemultcc'
        optimalvalue=retrievediffmodelinformation(ocStruct,'optimalvalue',field2arcidentifier(arcfield));
        inequalitycontrolconstraintnum=retrievediffmodelinformation(ocStruct,'inequalitycontrolconstraintnum');
        lagrangemultipliercontrolnametp1=retrievediffmodelinformation(ocStruct,'lagrangemultipliercontrolnametp1');
        for ii=1:inequalitycontrolconstraintnum.value
            if ii==1
                if inequalitycontrolconstraintnum.value>1
                    algebraic.term{ii}=['[' optimalvalue.value.(lagrangemultipliercontrolnametp1.value{ii}) '; ...'];
                else
                    algebraic.term{ii}=optimalvalue.value.(lagrangemultipliercontrolnametp1.value{ii});
                end
            elseif ii<inequalitycontrolconstraintnum.value
                algebraic.term{ii}=[optimalvalue.value.(lagrangemultipliercontrolnametp1.value{ii}) '; ...'];
            else
                algebraic.term{ii}=optimalvalue.value.(lagrangemultipliercontrolnametp1.value{ii});
            end
        end
        if inequalitycontrolconstraintnum.value>1
            algebraic.term{ii}=[algebraic.term{ii} ']'];
        end
    case 'constraint'
        inequalitycontrolconstraint=retrievediffmodelinformation(ocStruct,'inequalitycontrolconstraint');
        algebraic.term=inequalitycontrolconstraint.value;
    case 'symboliccontrolvalue'
        optimalvalue=retrievediffmodelinformation(ocStruct,'optimalvalue',field2arcidentifier(arcfield));
        controlnum=retrievediffmodelinformation(ocStruct,'controlnum');
        controlnametp1=retrievediffmodelinformation(ocStruct,'controlnametp1');
        for ii=1:controlnum.value
            if ii==1
                if controlnum.value>1
                    algebraic.term{ii}=['sym([''[' optimalvalue.value.(controlnametp1.value{ii}) ','' ...'];
                else
                    algebraic.term{ii}=['sym(''' optimalvalue.value.(controlnametp1.value{ii}) ''')'];
                end
            elseif ii<controlnum.value
                algebraic.term{ii}=['''' optimalvalue.value.(controlnametp1.value{ii}) ','' ...'];
            else
                algebraic.term{ii}=['''' optimalvalue.value.(controlnametp1.value{ii}) ']''])'];
            end
        end
    case 'symboliclagrangemultcc'
        optimalvalue=retrievediffmodelinformation(ocStruct,'optimalvalue',field2arcidentifier(arcfield));
        inequalitycontrolconstraintnum=retrievediffmodelinformation(ocStruct,'inequalitycontrolconstraintnum');
        lagrangemultipliercontrolnametp1=retrievediffmodelinformation(ocStruct,'lagrangemultipliercontrolnametp1');
        for ii=1:inequalitycontrolconstraintnum.value
            if ii==1
                if inequalitycontrolconstraintnum.value>1
                    algebraic.term{ii}=['sym([''[' optimalvalue.value.(lagrangemultipliercontrolnametp1.value{ii}) ','' ...'];
                else
                    algebraic.term{ii}=['sym(''' optimalvalue.value.(lagrangemultipliercontrolnametp1.value{ii}) ''')'];
                end
                
            elseif ii<inequalitycontrolconstraintnum.value
                algebraic.term{ii}=['''' optimalvalue.value.(lagrangemultipliercontrolnametp1.value{ii}) ','' ...'];
            else
                algebraic.term{ii}=['''' optimalvalue.value.(lagrangemultipliercontrolnametp1.value{ii}) ']''])'];
            end
        end
    case 'symboliccanonicalsystemmap'
        canonicalsystemimplicit=retrievediffmodelinformation(ocStruct,'canonicalsystemimplicit','','');
        for ii=1:numel(canonicalsystemimplicit.value)
            if ii==1
                algebraic.term{ii}=['sym([''[' canonicalsystemimplicit.value{ii} ';'' ...'];
            elseif ii<numel(canonicalsystemimplicit.value)
                algebraic.term{ii}=['''' canonicalsystemimplicit.value{ii} ';'' ...'];
            else
                algebraic.term{ii}=['''' canonicalsystemimplicit.value{ii} ']''])'];
            end
        end

    case 'symboliccanonicalsystemmapjacobian'
        arcidentifier=field2arcidentifier(arcfield);
        canonicalsystemjacobian=retrievediffmodelinformation(ocStruct,'canonicalsystemjacobian',arcidentifier,getsymkernel);
        numterm=numel(canonicalsystemjacobian.value);
        for ii=1:numterm
            if ii==1
                algebraic.term{ii}=['sym([''' canonicalsystemjacobian.value{ii} ';'' ...'];
            elseif ii<numterm
                algebraic.term{ii}=['''' canonicalsystemjacobian.value{ii} ';'' ...'];
            else
                algebraic.term{ii}=['''' canonicalsystemjacobian.value{ii} ''''];
            end
        end
        algebraic.term{ii}=[algebraic.term{ii} '])'];
    case 'adjointstateexplicit'
        arcidentifier=field2arcidentifier(arcfield);
        canonicalsystemexplicit=retrievediffmodelinformation(ocStruct,'adjointstateexplicit',arcidentifier,getsymkernel);
        numterm=numel(canonicalsystemexplicit.value);
        for ii=1:numterm
            if ii==1
                if numterm>1
                    algebraic.term{ii}=['[' canonicalsystemexplicit.value{ii} '; ...'];
                else
                    algebraic.term{ii}=canonicalsystemexplicit.value{ii};
                end
            elseif ii<numel(canonicalsystemexplicit.value)
                algebraic.term{ii}=[canonicalsystemexplicit.value{ii} '; ...'];
            else
                if numterm>1
                    algebraic.term{ii}=[canonicalsystemexplicit.value{ii} ']'];
                end
            end
        end

    case 'canonicalsystem'
        canonicalsystemimplicit=retrievediffmodelinformation(ocStruct,'canonicalsystemimplicit','',getsymkernel);
        numterm=numel(canonicalsystemimplicit.value);
        for ii=1:numterm
            if ii==1
                if numterm>1
                    algebraic.term{ii}=['[' canonicalsystemimplicit.value{ii} '; ...'];
                else
                    algebraic.term{ii}=canonicalsystemimplicit.value{ii};
                end
            elseif ii<numel(canonicalsystemimplicit.value)
                algebraic.term{ii}=[canonicalsystemimplicit.value{ii} '; ...'];
            else
                if numterm>1
                    algebraic.term{ii}=[canonicalsystemimplicit.value{ii} ']'];
                end
            end
        end

    case 'canonicalsystemjacobian'
        arcidentifier=field2arcidentifier(arcfield);
        canonicalsystemexplicitjacobian=retrievediffmodelinformation(ocStruct,'canonicalsystemjacobian',arcidentifier,getsymkernel);
        numterm=numel(canonicalsystemexplicitjacobian.value);
        for ii=1:numterm
            if ii==1
                if numterm>1
                    algebraic.term{ii}=['[' canonicalsystemexplicitjacobian.value{ii} '; ...'];
                else
                    algebraic.term{ii}=canonicalsystemexplicitjacobian.value{ii};
                end
            elseif ii<numel(canonicalsystemexplicitjacobian.value)
                algebraic.term{ii}=[canonicalsystemexplicitjacobian.value{ii} '; ...'];
            else
                if numterm>1
                    algebraic.term{ii}=[canonicalsystemexplicitjacobian.value{ii} ']'];
                end
            end
        end
    case 'canonicalsystemparameterjacobian'
        arcidentifier=field2arcidentifier(arcfield);
        canonicalsystemparameterjacobian=retrievediffmodelinformation(ocStruct,'canonicalsystemparameterjacobian',arcidentifier,getsymkernel);
        numterm=numel(canonicalsystemparameterjacobian.value);
        for ii=1:numterm
            if ii==1
                algebraic.term{ii}=['[' canonicalsystemparameterjacobian.value{ii} '; ...'];
            elseif ii<numterm
                algebraic.term{ii}=[canonicalsystemparameterjacobian.value{ii} '; ...'];
            else
                algebraic.term{ii}=canonicalsystemparameterjacobian.value{ii};
            end
        end
        algebraic.term{ii}=[algebraic.term{ii} ']'];
    case 'statedynamics'
        statedynamics=retrievediffmodelinformation(ocStruct,'statedynamics','','');
        numterm=numel(statedynamics.value);
        for ii=1:numel(statedynamics.value)
            if ii==1
                if numterm>1
                    algebraic.term{ii}=['[' statedynamics.value{ii} '; ...'];
                else
                    algebraic.term{ii}=statedynamics.value{ii};
                end
            elseif ii<numel(statedynamics.value)
                algebraic.term{ii}=[statedynamics.value{ii} '; ...'];
            else
                if numterm>1
                    algebraic.term{ii}=[statedynamics.value{ii} ']'];
                end
            end
        end
    case 'specificstatedynamics'
        arcidentifier=field2arcidentifier(arcfield);
        statedynamics=retrievediffmodelinformation(ocStruct,'specificstatedynamics',arcidentifier,getsymkernel);
        numterm=numel(statedynamics.value);
        for ii=1:numel(statedynamics.value)
            if ii==1
                if numterm>1
                    algebraic.term{ii}=['[' statedynamics.value{ii} '; ...'];
                else
                    algebraic.term{ii}=statedynamics.value{ii};
                end
            elseif ii<numel(statedynamics.value)
                algebraic.term{ii}=[statedynamics.value{ii} '; ...'];
            else
                if numterm>1
                    algebraic.term{ii}=[statedynamics.value{ii} ']'];
                end
            end
        end

    case 'specificstatedynamicsjacobian'
        arcidentifier=field2arcidentifier(arcfield);
        statedynamics=retrievediffmodelinformation(ocStruct,'specificstatedynamicsjacobian',arcidentifier,getsymkernel);
        numterm=numel(statedynamics.value);
        for ii=1:numel(statedynamics.value)
            if ii==1
                if numterm>1
                    algebraic.term{ii}=['[' statedynamics.value{ii} '; ...'];
                else
                    algebraic.term{ii}=statedynamics.value{ii};
                end
            elseif ii<numel(statedynamics.value)
                algebraic.term{ii}=[statedynamics.value{ii} '; ...'];
            else
                if numterm>1
                    algebraic.term{ii}=[statedynamics.value{ii} ']'];
                end
            end
        end

    case 'discobjectivefunctionjacobianx'
        arcidentifier=field2arcidentifier(arcfield);
        discobjectivefunctionjacobian=retrievediffmodelinformation(ocStruct,'discobjectivefunctionDX',arcidentifier,getsymkernel);
        statenum=retrievediffmodelinformation(ocStruct,'statenum');
        statecoord=1:statenum.value;
        for ii=statecoord
            if ii==1
                if statenum.value>1
                    algebraic.term{ii}=['[' discobjectivefunctionjacobian.value{ii} '; ...'];
                else
                    algebraic.term{ii}=[discobjectivefunctionjacobian.value{ii}];
                end
            elseif ii<statenum.value
                algebraic.term{ii}=[discobjectivefunctionjacobian.value{ii} '; ...'];
            else
                algebraic.term{ii}=discobjectivefunctionjacobian.value{ii};
            end
        end
        if statenum.value>1
            algebraic.term{ii}=[algebraic.term{ii} ']'];
        end
    case 'discobjectivefunctionjacobianl'
        arcidentifier=field2arcidentifier(arcfield);
        discobjectivefunctionjacobian=retrievediffmodelinformation(ocStruct,'discobjectivefunctionDX',arcidentifier,getsymkernel);
        statenum=retrievediffmodelinformation(ocStruct,'statenum');
        costatecoord=statenum.value+(1:statenum.value);
        counter=0;
        for ii=costatecoord
            counter=counter+1;
            if counter==1
                if statenum.value>1
                    algebraic.term{counter}=['[' discobjectivefunctionjacobian.value{ii} '; ...'];
                else
                    algebraic.term{counter}=[discobjectivefunctionjacobian.value{ii}];
                end
            elseif counter<statenum.value
                algebraic.term{counter}=[discobjectivefunctionjacobian.value{ii} '; ...'];
            else
                algebraic.term{counter}=discobjectivefunctionjacobian.value{ii};
            end
        end
        if statenum.value>1
            algebraic.term{counter}=[algebraic.term{counter} ']'];
        end

    case 'discobjectivefunctionparameterjacobian'
        arcidentifier=field2arcidentifier(arcfield);
        discobjectivefunctionjacobian=retrievediffmodelinformation(ocStruct,'discobjectivefunctionDP',arcidentifier,getsymkernel);
        numterm=numel(discobjectivefunctionjacobian.value);
        for ii=1:numterm
            if ii==1
                algebraic.term{ii}=['[' discobjectivefunctionjacobian.value{ii} ', ...'];
            elseif ii<numterm
                algebraic.term{ii}=[discobjectivefunctionjacobian.value{ii} ', ...'];
            else
                algebraic.term{ii}=discobjectivefunctionjacobian.value{ii};
            end
        end
        algebraic.term{ii}=[algebraic.term{ii} ']'];
    case 'discobjectivefunction'
        discobjectivefunction=retrievediffmodelinformation(ocStruct,'discobjectivefunction');
        algebraic.term{1}=discobjectivefunction.value;

    case 'parametervariables'
        parametername=retrieveodemodelinformation(ocStruct,'parametername');
        algebraic.term=parametername.value;

    case 'symbolicpontryaginfunction'
        pontryaginfunction=retrievediffmodelinformation(ocStruct,'pontryaginfunction');
        algebraic.term{1}=['sym(''' pontryaginfunction.value ''')'];

    case 'pontryaginfunction'
        pontryaginfunction=retrievediffmodelinformation(ocStruct,'pontryaginfunction');
        algebraic.term{1}=pontryaginfunction.value;
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