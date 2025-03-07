function arcinfo=cartel2dsymArcInfo(arcid)
% contains general arc information of the actual model
%
% this file was automatically created: 24-Aug-2024 15:01:54
% written by Dieter Grass, 2001 - 2024
	
switch arcid
	case 0
		arcinfo.ArcIdentifier=0; % unique identifier for the specific arc
	
		% information about the canonical system 
		arcinfo.statecoord=1:2;
		arcinfo.costatecoord=3:4;
		arcinfo.implicitcontrolcoord=5; % coordinate of the implicitly determined control values
		arcinfo.implicitcontrolindex=3; % index of the implicitly determined control values
	
		% the following information is redundant and has to be consistent
		% with the above definitions
		arcinfo.odedim=4; % dimension of the odes 
		arcinfo.aedim=1; % dimension of the algebraic equations
		arcinfo.daeorder=[1  1  1  1  0]; % order of the DAEs
	
	case 1
		arcinfo.ArcIdentifier=1; % unique identifier for the specific arc
	
		% information about the canonical system 
		arcinfo.statecoord=1:2;
		arcinfo.costatecoord=3:4;
		arcinfo.implicitcontrolcoord=5; % coordinate of the implicitly determined control values
		arcinfo.implicitcontrolindex=3; % index of the implicitly determined control values
	
		% the following information is redundant and has to be consistent
		% with the above definitions
		arcinfo.odedim=4; % dimension of the odes 
		arcinfo.aedim=1; % dimension of the algebraic equations
		arcinfo.daeorder=[1  1  1  1  0]; % order of the DAEs
	
	case 2
		arcinfo.ArcIdentifier=2; % unique identifier for the specific arc
	
		% information about the canonical system 
		arcinfo.statecoord=1:2;
		arcinfo.costatecoord=3:4;
		arcinfo.implicitcontrolcoord=5; % coordinate of the implicitly determined control values
		arcinfo.implicitcontrolindex=3; % index of the implicitly determined control values
	
		% the following information is redundant and has to be consistent
		% with the above definitions
		arcinfo.odedim=4; % dimension of the odes 
		arcinfo.aedim=1; % dimension of the algebraic equations
		arcinfo.daeorder=[1  1  1  1  0]; % order of the DAEs
	
	case 3
		arcinfo.ArcIdentifier=3; % unique identifier for the specific arc
	
		% information about the canonical system 
		arcinfo.statecoord=1:2;
		arcinfo.costatecoord=3:4;
		arcinfo.implicitcontrolcoord=[]; % coordinate of the implicitly determined control values
		arcinfo.implicitcontrolindex=[]; % index of the implicitly determined control values
	
		% the following information is redundant and has to be consistent
		% with the above definitions
		arcinfo.odedim=4; % dimension of the odes 
		arcinfo.aedim=0; % dimension of the algebraic equations
		arcinfo.daeorder=[1  1  1  1]; % order of the DAEs
	
	case 4
		arcinfo.ArcIdentifier=4; % unique identifier for the specific arc
	
		% information about the canonical system 
		arcinfo.statecoord=1:2;
		arcinfo.costatecoord=3:4;
		arcinfo.implicitcontrolcoord=[5  6]; % coordinate of the implicitly determined control values
		arcinfo.implicitcontrolindex=[1  2]; % index of the implicitly determined control values
	
		% the following information is redundant and has to be consistent
		% with the above definitions
		arcinfo.odedim=4; % dimension of the odes 
		arcinfo.aedim=2; % dimension of the algebraic equations
		arcinfo.daeorder=[1  1  1  1  0  0]; % order of the DAEs
	
	case 5
		arcinfo.ArcIdentifier=5; % unique identifier for the specific arc
	
		% information about the canonical system 
		arcinfo.statecoord=1:2;
		arcinfo.costatecoord=3:4;
		arcinfo.implicitcontrolcoord=5; % coordinate of the implicitly determined control values
		arcinfo.implicitcontrolindex=3; % index of the implicitly determined control values
	
		% the following information is redundant and has to be consistent
		% with the above definitions
		arcinfo.odedim=4; % dimension of the odes 
		arcinfo.aedim=1; % dimension of the algebraic equations
		arcinfo.daeorder=[1  1  1  1  0]; % order of the DAEs
	
	case 6
		arcinfo.ArcIdentifier=6; % unique identifier for the specific arc
	
		% information about the canonical system 
		arcinfo.statecoord=1:2;
		arcinfo.costatecoord=3:4;
		arcinfo.implicitcontrolcoord=[]; % coordinate of the implicitly determined control values
		arcinfo.implicitcontrolindex=[]; % index of the implicitly determined control values
	
		% the following information is redundant and has to be consistent
		% with the above definitions
		arcinfo.odedim=4; % dimension of the odes 
		arcinfo.aedim=0; % dimension of the algebraic equations
		arcinfo.daeorder=[1  1  1  1]; % order of the DAEs
	
	case 7
		arcinfo.ArcIdentifier=7; % unique identifier for the specific arc
	
		% information about the canonical system 
		arcinfo.statecoord=1:2;
		arcinfo.costatecoord=3:4;
		arcinfo.implicitcontrolcoord=5; % coordinate of the implicitly determined control values
		arcinfo.implicitcontrolindex=2; % index of the implicitly determined control values
	
		% the following information is redundant and has to be consistent
		% with the above definitions
		arcinfo.odedim=4; % dimension of the odes 
		arcinfo.aedim=1; % dimension of the algebraic equations
		arcinfo.daeorder=[1  1  1  1  0]; % order of the DAEs
	
	case 8
		arcinfo.ArcIdentifier=8; % unique identifier for the specific arc
	
		% information about the canonical system 
		arcinfo.statecoord=1:2;
		arcinfo.costatecoord=3:4;
		arcinfo.implicitcontrolcoord=[]; % coordinate of the implicitly determined control values
		arcinfo.implicitcontrolindex=[]; % index of the implicitly determined control values
	
		% the following information is redundant and has to be consistent
		% with the above definitions
		arcinfo.odedim=4; % dimension of the odes 
		arcinfo.aedim=0; % dimension of the algebraic equations
		arcinfo.daeorder=[1  1  1  1]; % order of the DAEs
	
	case 9
		arcinfo.ArcIdentifier=9; % unique identifier for the specific arc
	
		% information about the canonical system 
		arcinfo.statecoord=1:2;
		arcinfo.costatecoord=3:4;
		arcinfo.implicitcontrolcoord=5; % coordinate of the implicitly determined control values
		arcinfo.implicitcontrolindex=1; % index of the implicitly determined control values
	
		% the following information is redundant and has to be consistent
		% with the above definitions
		arcinfo.odedim=4; % dimension of the odes 
		arcinfo.aedim=1; % dimension of the algebraic equations
		arcinfo.daeorder=[1  1  1  1  0]; % order of the DAEs
	
	case 10
		arcinfo.ArcIdentifier=10; % unique identifier for the specific arc
	
		% information about the canonical system 
		arcinfo.statecoord=1:2;
		arcinfo.costatecoord=3:4;
		arcinfo.implicitcontrolcoord=5; % coordinate of the implicitly determined control values
		arcinfo.implicitcontrolindex=2; % index of the implicitly determined control values
	
		% the following information is redundant and has to be consistent
		% with the above definitions
		arcinfo.odedim=4; % dimension of the odes 
		arcinfo.aedim=1; % dimension of the algebraic equations
		arcinfo.daeorder=[1  1  1  1  0]; % order of the DAEs
	
	case 11
		arcinfo.ArcIdentifier=11; % unique identifier for the specific arc
	
		% information about the canonical system 
		arcinfo.statecoord=1:2;
		arcinfo.costatecoord=3:4;
		arcinfo.implicitcontrolcoord=[]; % coordinate of the implicitly determined control values
		arcinfo.implicitcontrolindex=[]; % index of the implicitly determined control values
	
		% the following information is redundant and has to be consistent
		% with the above definitions
		arcinfo.odedim=4; % dimension of the odes 
		arcinfo.aedim=0; % dimension of the algebraic equations
		arcinfo.daeorder=[1  1  1  1]; % order of the DAEs
	
	case 12
		arcinfo.ArcIdentifier=12; % unique identifier for the specific arc
	
		% information about the canonical system 
		arcinfo.statecoord=1:2;
		arcinfo.costatecoord=3:4;
		arcinfo.implicitcontrolcoord=[]; % coordinate of the implicitly determined control values
		arcinfo.implicitcontrolindex=[]; % index of the implicitly determined control values
	
		% the following information is redundant and has to be consistent
		% with the above definitions
		arcinfo.odedim=4; % dimension of the odes 
		arcinfo.aedim=0; % dimension of the algebraic equations
		arcinfo.daeorder=[1  1  1  1]; % order of the DAEs
	
	case 13
		arcinfo.ArcIdentifier=13; % unique identifier for the specific arc
	
		% information about the canonical system 
		arcinfo.statecoord=1:2;
		arcinfo.costatecoord=3:4;
		arcinfo.implicitcontrolcoord=[]; % coordinate of the implicitly determined control values
		arcinfo.implicitcontrolindex=[]; % index of the implicitly determined control values
	
		% the following information is redundant and has to be consistent
		% with the above definitions
		arcinfo.odedim=4; % dimension of the odes 
		arcinfo.aedim=0; % dimension of the algebraic equations
		arcinfo.daeorder=[1  1  1  1]; % order of the DAEs
	
	case 14
		arcinfo.ArcIdentifier=14; % unique identifier for the specific arc
	
		% information about the canonical system 
		arcinfo.statecoord=1:2;
		arcinfo.costatecoord=3:4;
		arcinfo.implicitcontrolcoord=[]; % coordinate of the implicitly determined control values
		arcinfo.implicitcontrolindex=[]; % index of the implicitly determined control values
	
		% the following information is redundant and has to be consistent
		% with the above definitions
		arcinfo.odedim=4; % dimension of the odes 
		arcinfo.aedim=0; % dimension of the algebraic equations
		arcinfo.daeorder=[1  1  1  1]; % order of the DAEs
	
end
