function out = harvest2DHybridStructureDesign
% the names guard and reset are taken from the hybrid literature as it can
% be found, e.g., in A.J. van der Schaft and H. Schumacher (1999), An
% Introduction to Hybrid Dynamical Systems, Springer-Verlag

out{1}=@hybridinfo;
out{2}=@domain;
out{3}=@hybriddynamics;
out{4}=@guard;
out{5}=@reset;
out{6}=@switchtime;
out{7}=@hybridjacobian;
out{8}=[];%@guardjacobian;
out{9}=@resetjacobian;

%--------------------------------------------------------------------------
function infoH=hybridinfo()

infoH.statuslabel=[1:2]; % replaces name of arcnum arcidentifier
infoH.edge={'1 2','2 1'}; % a cell of strings which describes the possible transitions between the states of the hybrid system    

function infoS=domain(statuslabel)
% returns info about the domain and its used discretization for each of the hybrid statuses
switch statuslabel
    case 1
        infoS.statuslabel=1; % unique identifier for the specific arc

        % information about the canonical system 
        infoS.statecoord=1:2; % coordinate of the states in the system of DAEs or ODEs, respectively 
        infoS.costatecoord=3:4; % coordinate of the costates in the system of DAEs or ODEs, respectively
        infoS.implicitcontrolcoord=5; % coordinate of the implicitly determined control values
        infoS.implicitcontrolindex=1; % index of the implicitly determined control values

        % the following information is redundant and has to be consistent
        % with the above definitions
        infoS.odedim=4; % dimension of the odes 
        infoS.aedim=1; % dimension of the algebraic equations
        infoS.odeorder=[1 1 1 1 0]; % order of the DAEs
        
        % general information for the discretization of each arc (in case
        % that the SBVPOC Solver is used)
        infoS.discretization.collocationmethod='g'; % possible collocation methods are 'Gaussian' or 'Lobatto'
        infoS.discretization.numcollocationpoints=6; % number of collocation points
    case 2
        infoS.statuslabel=2; % unique identifier for the specific arc

        % information about the canonical system 
        infoS.statecoord=1:2; % coordinate of the states in the system of DAEs or ODEs, respectively
        infoS.costatecoord=3:4; % coordinate of the costates in the system of DAEs or ODEs, respectively
        infoS.implicitcontrolcoord=[]; % coordinate of the implicitly determined control values
        infoS.implicitcontrolindex=[]; % index of the implicitly determined control values

        % the following information is redundant and has to be consistent
        % with the above definitions
        infoS.odedim=4; % dimension of the odes 
        infoS.aedim=[]; % dimension of the algebraic equations
        infoS.odeorder=[1 1 1 1]; % order of the DAEs
        
        % general information for the discretization of each arc (in case
        % that the SBVPOC Solver is used)
        infoS.discretization.collocationmethod='g'; % possible collocation methods are 'Gaussian' or 'Lobatto'
        infoS.discretization.numcollocationpoints=4; % number of collocation points
end

%--------------------------------------------------------------------------
function dxdt=hybriddynamics(t,dynVar,pararg,statusidx)

dxdt=harvest2DCanonicalSystem(t,dynVar,pararg,statusidx);

%--------------------------------------------------------------------------
function val=guard(t,dynVar,pararg,edge)
% condition on the (continuous) variables (states, costates) when a
% switching from status1 to status2 takes place. edge='status 1 status2'
% this condition determines the switching time
switch edge
    case '1 2'
		mc=harvest2DConstraint([],dynVar,pararg,1);
		val=mc(1);
    case '2 1'
        lmmc=harvest2DLagrangeMultiplier([],dynVar,pararg,1);
		val=lmmc(1);
end

%--------------------------------------------------------------------------
function val=reset(t,dynVarSi,dynVarSj,pararg,edge)
%transformation of the (continuous) variables (states, costates) when
%switching switching from status1 to status2. edge='status 1 status2' 

val=dynVarSi-dynVarSj;


%--------------------------------------------------------------------------
function val=switchtime(t,dynVar,pararg,edge)
% here a fixed switchtime for a specific edge can be specified, otherwsie
% it is determined by the guard
val=[];

%--------------------------------------------------------------------------
function [jac jacp]=hybridjacobian(t,dynVar,pararg,statusidx)

jac=harvest2DCanonicalSystemJacobian(t,dynVar,pararg,statusidx);
jacp=harvest2DParameterJacobian(t,dynVar,pararg,statusidx);

%--------------------------------------------------------------------------
function [jac jacp]=guardjacobian(t,dynVar,pararg,statusidx)
jac=[];
jacp=[];

%--------------------------------------------------------------------------
function [jaci jacj jacp]=resetjacobian(t,dynVarSi,dynVarSj,pararg,edge)

jaci=[1 0;0 1];
jacj=[-1 0;0 -1];
jacp=[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; ...
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
