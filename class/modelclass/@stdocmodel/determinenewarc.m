function newarcStruct=determinenewarc(ocObj,violationStruct)
%
% DETERMINENEWARC returns information about an arc at which a violation
% occurs. 
%
% DETERMINENEWARC(OCOBJ,VIOLATIONSTRUCT) it is assumed that only one
% violation occurs at a specific grid point. If multiple violations occur
% this usually refer to a numerical artefact and has to be processedby the
% user. The violation structure VIOLATIONSTRUCT has to consist of the following
% fields: 
% arcarg ... the arc argument for which the violation occurs
% arcnum ... the number of the arc at which the violation occurs
% violationmat ... a 0,1 matrix, where 0 corresponds to no violation and 1
%   that a violation occurs. The size of violationmat is given by two times
%   the number of (inequality) constraints (M) and the number of grid
%   points. Values 1 occuring at row j<=M denote a violation of the
%   constraint j, whereas M<j<=2*M denote that the jth constraint becomes
%   inactive. 
% constraintvalue ... is of the same size as violationmat and contains the
%   actual values of the constraint.
% rows and cols ... refer to the row and column for the violationmat where
%   the violation occurs 
%
% NEWARCSTRUCT=DETERMINENEWARC(OCOBJ,VIOLATIONSTRUCT) NEWARCSTRUCT consists
% of the follwoing fields:
% newarcarg ... the arc argument for the new arc taking the violation into
%   account
% violationposition ... each column consists of the left and right value of
%   the interval with consecutive values of cols e.g. 
%   cols=[2 3 4 5 8 9 11] => violationposition=[[2;5] [8;9] [11;11]]
% constraintindex ... for each violation interval the index of the
%   corresponding violation (rows value) is returned, e.g.
%   
%
% Example:
% violationStruct.cols[2 3 4 5 8 9 11];
% violationStruct.rows=[1 1 1 1 4 4 4 1];
% yields
% newarcstruct.violationposition=[[2;5] [8;9] [11;11]]
% newarcstruct.constraintindex=[1 4 1]
%
% the following examples will not be processed
%
% violationStruct.cols[2 2 4 5 8 9 11]; multiple violations at the same
% grid point
%
% violationStruct.cols[2 3 4 5 8 9 11];
% violationStruct.rows=[1 1 1 1 4 3 4 1]; different violations within one
% violation interval 

newarcStruct=[];
if isempty(violationStruct)
    return
end
if isempty(ocObj)
    return
end
if ~isstruct(violationStruct) || ~isfield(violationStruct,'arcarg') || ~isfield(violationStruct,'rows') ...
         || ~isfield(violationStruct,'cols') || ~isfield(violationStruct,'violationmat') ...
         || ~isfield(violationStruct,'constraintvalue') || ~isfield(violationStruct,'arcnum')
    ocmaterror('Second input argument is not a violation structure')
end
if numel(violationStruct.cols)~=numel(violationStruct.rows)
    ocmaterror('Number of violation columns and rows are different.')
end
if numel(violationStruct.cols)~=numel(unique(violationStruct.cols))
    % usually violations do not occur at the same point
    % therefore the generation of an initial solution is
    % left to the user.
    ocmatmsg('Multiple violations occur at the same grid point.\nThis case will not be processed automatically.\n')
    return
end
constraintnum=controlconstraintnum(ocObj)+stateconstraintnum(ocObj);
leftindex=[1 find(diff(violationStruct.cols)>1)+1];
rightindex=[find(diff(violationStruct.cols)>1) numel(violationStruct.cols)];
violationposition=[leftindex;rightindex];
numviolation=numel(leftindex);
newarcarg=zeros(1,numviolation);
constraintindex=zeros(1,numviolation);
for ii=1:numviolation
    % test if the same violation occurs along the actual violation interval,
    % a violation interval is defined as an interval with consecutive
    % values for the columns, e.g. cols=[2 3 4 5 8 9 11] consists of three
    % violation intervals [[2;5] [8;9] [11;11]]
    % rows than has to be of the structure rows=[idx1 idx1 idx1 idx1 idx2 idx2 idx3] 
    % which means that at each violation interval the same violation
    % occurs. Only this case will be processed automatically. In the other
    % case the user has to determine the new solution structure her/himself.
    if ~all(violationStruct.rows(leftindex(ii):rightindex(ii))==violationStruct.rows(leftindex(ii)))
        newarcStruct=[];
        ocmatmsg('Different violations occur at the same violation interval.\nThis case will not be processed automatically.\n')
        return
    end
    % idx ... is a vector of the length of control constraints, 0 ..
    % constraint inactive, 1 constraint active. This vector corresponds
    % ont-to-one to an arc argument.
    idx=arcarg2constraintcombinationindex(ocObj,violationStruct.arcarg);
    newidx=idx;
    if violationStruct.rows(leftindex(ii))<=constraintnum
        % a constraint is violated and becomes active, corresponding
        % position at newidx is set to 1
        newidx(violationStruct.rows(leftindex(ii)))=1;
    elseif violationStruct.rows(leftindex(ii))<=2*constraintnum
        % the Lagrange multiplier is negative and the constraint
        % becomes inactive, corresponding position at newidx is set to 0
        newidx(violationStruct.rows(leftindex(ii))-constraintnum)=0;
    end
    tmparg=constraintcombinationindex2arcarg(ocObj,newidx);
    if isempty(tmparg)
        newarcStruct=[];
        ocmatmsg('The constraint combination for index vector ''%d'' is unknown.\nThe user should add this case to the initalization file and rerun the initialization process.\n',newidx)
        return
    end
    newarcarg(ii)=tmparg;
    constraintindex(ii)=violationStruct.rows(leftindex(ii));
end
newarcStruct.newarcarg=newarcarg; % the new arc argument for the arc where the violation occurs
newarcStruct.violationposition=violationStruct.cols(violationposition); % the columns along the violation occurs
newarcStruct.constraintindex=constraintindex;