function arcInfo = harvest2DArcInfo(arcarg)
%
switch arcarg
    case 1
        arcInfo.ArcIdentifier=1; % unique identifier for the specific arc

        % information about the canonical system 
        arcInfo.StateCoord=1:2; % coordinate of the states in the system of DAEs or ODEs, respectively
        arcInfo.CostateCoord=3:4; % coordinate of the costates in the system of DAEs or ODEs, respectively
        arcInfo.ImplicitControlCoord=5; % coordinate of the implicitly determined control values
        arcInfo.ImplicitControlIndex=1; % index of the implicitly determined control values

        % the following information is redundant and has to be consistent
        % with the above definitions
        arcInfo.ODEDim=4; % dimension of the odes 
        arcInfo.AEDim=1; % dimension of the algebraic equations
        arcInfo.ODEOrder=[1 1 1 1 0]; % order of the DAEs
        
        % general information for the discretization of each arc (in case
        % that the SBVPOC Solver is used)
        arcInfo.CollocationMethod='g'; % possible collocation methods are 'Gaussian' or 'Lobatto'
        arcInfo.NumCollocationPoints=6; % number of collocation points
    case 2
        arcInfo.ArcIdentifier=2; % unique identifier for the specific arc

        % information about the canonical system 
        arcInfo.StateCoord=1:2; % coordinate of the states in the system of DAEs or ODEs, respectively
        arcInfo.CostateCoord=3:4; % coordinate of the costates in the system of DAEs or ODEs, respectively
        arcInfo.ImplicitControlCoord=[]; % coordinate of the implicitly determined control values
        arcInfo.ImplicitControlIndex=[]; % index of the implicitly determined control values

        % the following information is redundant and has to be consistent
        % with the above definitions
        arcInfo.ODEDim=4; % dimension of the odes 
        arcInfo.AEDim=[]; % dimension of the algebraic equations
        arcInfo.ODEOrder=[1 1 1 1]; % order of the DAEs
        
        % general information for the discretization of each arc (in case
        % that the SBVPOC Solver is used)
        arcInfo.CollocationMethod='g'; % possible collocation methods are 'Gaussian' or 'Lobatto'
        arcInfo.NumCollocationPoints=4; % number of collocation points
end
