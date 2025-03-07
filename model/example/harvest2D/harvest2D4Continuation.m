function out=harvest2D4Continuation()

out{1}=@canonicalsystem;
out{2}=@jacobian;
out{3}{1}=@bcinitial;
out{3}{2}=@bcasymptotic;
out{3}{3}=@bctransversality;
out{4}{1}=@bcjacobianinitial;
out{4}{2}=@bcjacobianasymptotic;
out{4}{3}=@bcjacobiantransversality;
% functions for hybrid structure and its discretization
out{5}{1}=@hybridinfo;
out{5}{2}=@domain;
out{5}{3}=@guard;
out{5}{4}=@reset;
out{5}{5}=@switchtime;
out{5}{7}=@jacobianguard;
out{5}{8}=@jacobianreset;
out{5}{9}=@domaindiscretization;
out{5}{10}=@timesettransformation;

out{11}=@plotcontinuation;
out{12}=@testadmissibility;

out{20}=@datapath;

function dxdt=canonicalsystem(t,dynVar,pararg,status)
%dxdt=harvest2DCanonicalSystem(t,dynVar,pararg,status);
dxdt=harvest2DCanonicalSystemIC(t,dynVar,pararg,status);

%--------------------------------------------------------------------------
function J=jacobian(t,dynVar,pararg,status)
J=harvest2DCanonicalSystemICJacobian(t,dynVar,pararg,status);

%--------------------------------------------------------------------------
function val=guard(dynVaral,dynVarar,dynVarbl,dynVarbr,pararg,freepar,status,edge)
% condition on the (continuous) variables (states, costates) when a
% switching from status1 to status2 takes place. edge=[status 1 status2]
% this condition determines the switching time
if isempty(edge)
    val=[];
else
    val=harvest2DHamiltonian([],dynVarbl,pararg,edge(1))-harvest2DHamiltonian([],dynVarbr,pararg,edge(2));
end

%--------------------------------------------------------------------------
function val=reset(dynVaral,dynVarar,dynVarbl,dynVarbr,pararg,freepar,status,edge)
%transformation of the (continuous) variables (states, costates) when
%switching switching from status1 to status2. edge=[status 1 status2] 

if isempty(edge)
    val=[];
else
    val=dynVaral(1:4)-dynVarar(1:4);
end

%--------------------------------------------------------------------------
function res=bcinitial(dynVaral,dynVarar,dynVarbl,dynVarbr,pararg,freepar,status,edge,targetcoordinate,continuationvector,initialstate)
res=dynVarar(targetcoordinate)-initialstate;

%--------------------------------------------------------------------------
function res=bcasymptotic(dynVaral,dynVarar,dynVarbl,dynVarbr,pararg,freepar,status,edge,asymptoticmatrix,saddlepoint)
res=asymptoticmatrix*(dynVarbl(1:4)-saddlepoint);

%--------------------------------------------------------------------------
function [Jal Jar Jbl Jbr Jpar]=bcjacobianinitial(dynVaral,dynVarar,dynVarbl,dynVarbr,pararg,freepar,status,edge,targetcoordinate,continuationvector,initialstate)
Jal=[];
Jbr=[];
switch status
    case 1
        if numel(targetcoordinate)==2
            Jar=[0 0 0 0 0;0 0 0 0 0];
            Jar(targetcoordinate(1),targetcoordinate(1))=1;
            Jar(targetcoordinate(2),targetcoordinate(2))=1;
            Jbl=[0 0 0 0 0;0 0 0 0 0];
            if ~isempty(edge)
                Jbr=[0 0 0 0;0 0 0 0];
            end
        else
        end
    case 2
        if numel(targetcoordinate)==2
            Jar=[0 0 0 0;0 0 0 0];
            Jar(targetcoordinate(1),targetcoordinate(1))=1;
            Jar(targetcoordinate(2),targetcoordinate(2))=1;
            Jbl=[0 0 0 0;0 0 0 0];
            if ~isempty(edge)
                Jbr=[0 0 0 0 0;0 0 0 0 0];
            end
        else
        end
end
if numel(targetcoordinate)==2
    numfree=numel(freepar);
    Jpar=zeros(2,numfree);
    Jpar(1:2,numfree)=-continuationvector;
end
%--------------------------------------------------------------------------
function [Jal Jar Jbl Jbr Jpar]=bcjacobianasymptotic(dynVaral,dynVarar,dynVarbl,dynVarbr,pararg,freepar,status,edge,asymptoticmatrix,saddlepoint)
Jal=[];
Jbr=[];
switch status
    case 1
        if size(asymptoticmatrix,1)==2
            Jar=[0 0 0 0 0;0 0 0 0 0];
            Jbl=[0 0 0 0 0;0 0 0 0 0];
            Jbl(:,1:4)=asymptoticmatrix;
            Jpar=zeros(2,numel(freepar));
            if ~isempty(edge)
                Jal=[0 0 0 0;0 0 0 0];
            end
        else
        end
    case 2
        if size(asymptoticmatrix,1)==2
            Jar=[0 0 0 0;0 0 0 0];
            Jbl=asymptoticmatrix;
            Jpar=zeros(2,numel(freepar));
            if ~isempty(edge)
                Jal=[0 0 0 0 0;0 0 0 0 0];
            end
        else
        end
end

%--------------------------------------------------------------------------
function [Jal Jar Jbl Jbr Jpar]=jacobianguard(dynVaral,dynVarar,dynVarbl,dynVarbr,pararg,freepar,status,edge)
Jar=[];
Jal=[];
Jbr=[];
Jbl=[];
Jpar=zeros(1,numel(freepar));
switch edge(1)
    case 1
        if ~isempty(dynVaral)
            Jal=[0 0 0 0];
        end
        Jar=[0 0 0 0 0];
        Jbl=harvest2DDHamiltonianDx([],dynVarbl,pararg,edge(1));
        Jbr=-harvest2DDHamiltonianDx([],dynVarbr,pararg,edge(2));
    case 2
        if ~isempty(dynVaral)
            Jal=[0 0 0 0 0];
        end
        Jar=[0 0 0 0];
        Jbl=harvest2DDHamiltonianDx([],dynVarbl,pararg,edge(1));
        Jbr=-harvest2DDHamiltonianDx([],dynVarbr,pararg,edge(2));
end

%--------------------------------------------------------------------------
function [Jal Jar Jbl Jbr Jpar]=jacobianreset(dynVaral,dynVarar,dynVarbl,dynVarbr,pararg,freepar,status,edge)
Jbr=[];
Jpar=zeros(4,numel(freepar));
switch edge(1)
    case 1
        eye4=eye(4);
        Jal=zeros(5);
        Jal(1:4,1:4)=eye4;
        Jar=-eye4;
        Jbl=zeros(4);
        if ~isempty(dynVarbr)
            Jbr=zeros(4,5);
        end
    case 2
        eye4=eye(4);
        Jar=zeros(5);
        Jar(1:4,1:4)=-eye4;
        Jal=eye4;
        Jbl=zeros(4,5);
        if ~isempty(dynVarbr)
            Jbr=zeros(4);
        end
end
%--------------------------------------------------------------------------
function infoH=hybridinfo()

infoH.statuslabel=[1:2]; % replaces name of arcnum arcidentifier
infoH.edge=[1 2;2 1]; % each row represents a possible transition between the states of the hybrid system    

%--------------------------------------------------------------------------
function infoD=domain(statuslabel)
% returns info about the domain and its used discretization for each of the hybrid statuses
switch statuslabel
    case 1
        infoD.statuslabel=1; % unique identifier for the specific arc

        % information about the canonical system 
        infoD.statecoord=1:2; % coordinate of the states in the system of DAEs or ODEs, respectively 
        infoD.costatecoord=3:4; % coordinate of the costates in the system of DAEs or ODEs, respectively
        infoD.implicitcontrolcoord=5; % coordinate of the implicitly determined control values
        infoD.implicitcontrolindex=1; % index of the implicitly determined control values

        % the following information is redundant and has to be consistent
        % with the above definitions
        infoD.odedim=4; % dimension of the odes 
        infoD.aedim=1; % dimension of the algebraic equations
        infoD.odeorder=[1 1 1 1 0]; % order of the DAEs
        
    case 2
        infoD.statuslabel=2; % unique identifier for the specific arc

        % information about the canonical system 
        infoD.statecoord=1:2; % coordinate of the states in the system of DAEs or ODEs, respectively
        infoD.costatecoord=3:4; % coordinate of the costates in the system of DAEs or ODEs, respectively
        infoD.implicitcontrolcoord=[]; % coordinate of the implicitly determined control values
        infoD.implicitcontrolindex=[]; % index of the implicitly determined control values

        % the following information is redundant and has to be consistent
        % with the above definitions
        infoD.odedim=4; % dimension of the odes 
        infoD.aedim=[]; % dimension of the algebraic equations
        infoD.odeorder=[1 1 1 1]; % order of the DAEs
end

%--------------------------------------------------------------------------
function val=switchtime(edge)
% here a fixed switchtime for a specific edge can be specified, otherwsie
% it is determined by the guard
val=[];

%--------------------------------------------------------------------------
function infoD=domaindiscretization(statuslabel)
% returns info about the domain and its used discretization for each of the hybrid statuses
switch statuslabel
    case 1
        infoD.collocationmethod='g'; % possible collocation methods are 'Gaussian' or 'Lobatto'
        infoD.numcollocationpoints=4; % number of collocation points
        infoD.equilibrium.initnummesh=40;
    case 2
        infoD.collocationmethod='g'; % possible collocation methods are 'Gaussian' or 'Lobatto'
        infoD.numcollocationpoints=2; % number of collocation points
        infoD.equilibrium.initnummesh=40;
end

%--------------------------------------------------------------------------
function infoT=timesettransformation(statuslabel)
infoT.normalization=1; % the nth timenterval is normalized to [n,n+1]
infoT.infinity=-0.1; % for an infinite time horizon problem the last interval [N,inf)->[N,N+1]
infoT.asymptoticapproximation=1;
% infoT.infinity=0; % for an infinite time horizon problem the last interval [N,inf)->[N,N+1], infoT.infinity<0
% infoT.asymptoticapproximation=25;

%--------------------------------------------------------------------------
function h=plotcontinuation(t,dynVar,modelpar,arcarg,freepar,tangent,plotflag)
h=harvest2DPlotContinuation(t,dynVar,modelpar,arcarg,freepar,tangent,plotflag);

function [out labelS]=testadmissibility(t,dynVar,pararg,arcarg)
labelS(1).aracarg=1;
labelS(1).info={'mc1','lm1'};
labelS(2).aracarg=2;
labelS(2).info={'mc1','lm1'};

out=[harvest2DConstraint(t,dynVar,pararg,arcarg); ...
    harvest2DLagrangeMultiplier(t,dynVar,pararg,arcarg)]+1e-4;

%--------------------------------------------------------------------------
function pathname=datapath()

pathname='F:\home\dieter\Eigene Dateien\Matlab\toolbox\numtools\ocmatnew\model\harvest2D\data';
