function ocEP=calcic(ocObj,varargin)
%
% CALCIC calculates implicit control for given state-costate values.
%

depvar0=[]; % initial guess for equilibria
opt=[]; % option for fsolve
arcarg=[];
equationfile='';
jacobianfile='';
optimizationfile='';
gradientfile='';
method='';

if isempty(ocObj)
    return
end
if nargin>=2
    depvar0=varargin{1};
end
if nargin>=3
    arcarg=varargin{2};
end
if nargin>=4
    opt=varargin{3};
end
if isempty(opt)
    opt=defaultocoptions;
end
if isempty(arcarg)
    arcarg=arcargument(ocObj);
end
for ii=4:2:nargin-1
    if ischar(varargin{ii})
        eval([varargin{ii} '=varargin{ii+1};'])
    end
end
if isempty(method)
    method='zero';
end
switch method
    case 'zero'
        if isempty(equationfile)
            equationfile='AlgebraicEquationI';
        end
        if isempty(jacobianfile)
            jacobianfile='AlgebraicEquationJacobian';
        end
    case 'grad'
        if isempty(optimizationfile)
            optimizationfile='Hamiltonian';
        end
        if isempty(gradientfile)
            gradientfile='DHamiltonianDU';
        end
        
        [A,b,lb,ub]=generateconstraintvariable(ocObj);
end
MessageDisplay=strcmp(getocoptions(opt,'INIT','MessageDisplay'),'on');
numarc=numel(arcarg);
if MessageDisplay
    ocmatmsg('\nStart searching for implicit control values.\n')
end
counter=0;
statecostatecoord=1:2*statenum(ocObj);
numinitpt=size(depvar0,2);
numep=numarc*numinitpt;
epstruct=struct('y',cell(1,numep),'arcarg', cell(1,numep),'solverinfo', cell(1,numep));
epcounter=0;
while 1
    counter=counter+1;
    if counter>numarc
        break
    end
    act_arcarg=arcarg(counter);
    if MessageDisplay
        ocmatmsg('\nActual arc identifier for calculation %d\n',act_arcarg)
    end
    for ii=1:numinitpt
        coord=implicitcontrolcoordinate(ocObj,act_arcarg);
        coord=(1:length(coord))+statecostatecoord(end);
        epcounter=epcounter+1;
        switch method
            case 'zero'
                [ctrlval,fval,exitflag,info,jacob]=ocmateqsolve(ocObj,equationfile,0,depvar0(coord,ii),act_arcarg,opt,[],depvar0(statecostatecoord,ii));
            case 'grad'
                optimStruct.optimizationfile=optimizationfile;
                optimStruct.gradientfile=gradientfile;
                optimStruct.lowerbound=lb;
                optimStruct.upperbound=ub;
                optimStruct.mixedconstraintb=b;
                optimStruct.mixedconstraintA=A;
                optimStruct.option=opt;
                optimStruct.statecostate=depvar0(statecostatecoord,ii);
                optimStruct.time=0;
                optimStruct.arcarg=act_arcarg;
                [ctrlval,fval,exitflag,info,jacob]=ocmatoptimize(ocObj,depvar0(coord,ii),optimStruct);
        end

        epstruct(epcounter).y=[depvar0(statecostatecoord,ii);ctrlval];
        epstruct(epcounter).arcarg=act_arcarg;
        epstruct(epcounter).solverinfo=info;
        epstruct(epcounter).solverinfo.exitflag=exitflag;
        epstruct(epcounter).solverinfo.fval=fval;
        epstruct(epcounter).solverinfo.jacobian=jacob;
    end

end
epstruct(epcounter+1:numep)=[];

% calculate linearization and generate dynprimitive objects
ocEP=cell(epcounter,1);
for ii=1:epcounter
    ocEP{ii}=dynprimitive(epstruct(ii),ocObj);
    ocEP{ii}=gdynprimitive(ocEP{ii});
end

function [A,b,lb,ub]=generateconstraintvariable(ocObj,x)
A=[];
b=[];
lb=[];
ub=[];

if nargin==1
    x=[];
end
cstr=constraint(ocObj);
ctrlvariable=controlvariable(ocObj);
if isempty(cstr)
    return
end
cstrs=subsparametervalue(ocObj,cstr);
dcstrdu=subsparametervalue(ocObj,mystr2sym(removematrixstring(ocmatjacobian(removematrixstring(char(cstr.')),cell2vectorstring(ctrlvariable),getsymkernel))));

lb=repmat(-inf,length(ctrlvariable),1);
ub=repmat(inf,length(ctrlvariable),1);
A=zeros(length(cstr),length(ctrlvariable));
b=zeros(length(cstr),1);
for ii=1:size(dcstrdu,1)
    idx=find(dcstrdu(ii,:));
    if ~isempty(idx)
        if length(idx)==1
            val=subs(cstrs(ii),ctrlvariable{idx},0);
            if val<0
                ub(idx)=-val;
            else
                lb(idx)=val;
            end
        else
            ocmatmsg('Not implemented yet.\n')
%             A(ii,idx)=1;
%             for jj=idx
%                 if val<0
%                     b(idx)=-val;
%                 else
%                     b(idx)=val;
%                 end
%             end
        end
    end
end
if sum(abs(A(:)))==0
    A=[];
    b=[];
end
if all(isinf(lb))
    lb=[];
end
if all(isinf(ub))
    ub=[];
end

