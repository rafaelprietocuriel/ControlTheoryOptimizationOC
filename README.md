# The Mexican dilemma between social and security
Code created to run in Matlab to obtain the best possible policy strategies to reduce cartel violence and cartel harm.

#### software
Matlab R2007b and newer
tested with R2007b and R2023b

#### Instructions
% Place the 'ocmat' folder on your local system.
% Change to this folder from the MATLAB command window.

>> init  % Adds the necessary folders to the MATLAB path

% Set the options for the continuation process.
% This defines various parameters used for solving boundary value problems and continuation.
>> opt = setocoptions('OCCONTARG', 'MinStepWidth', 5e-2, 'MaxStepWidth', 5e5, ...
    'InitStepWidth', 1e-2, 'TotalRelativeDistance', 1, 'CheckAdmissibility', 1, ...
    'MaxDistance', 5e2, 'MaxContinuationSteps', inf, 'SBVPOC', 'MeshAdaptAbsTol', 1e-6, ...
    'MeshAdaptRelTol', 1e-4, 'MeshAdaptation', 1, 'NEWTON', 'MaxNewtonIters', 15, ...
    'MaxProbes', 10, 'CheckSingular', 0, 'RelTol', 1e-5, 'BVP', 'FJacobian', 1, ...
    'BCJacobian', 0, 'Vectorized', 'on', 'GENERAL', 'AdmissibleTolerance', 1e-4, ...
    'ZeroDeviationTolerance', 5e-3, 'BVPMethod', 'gbvp4c', 'INIT', 'TestModus', 'off');

>> m = stdocmodel('cartel2dsym0');  % Load the model named 'cartel2dsym0' from the data folder.

% Load the model with specific parameter values from the 'data' folder of the model
>> load(m, [], 1, 'BaseModel');

% Provide initial values for the equilibria computation. These values are from a previous numerical analysis.
>> X = [36472, 36472, -3679734, -3679734, 34400330, 34400330];

% Calculate the equilibrium using the given initial values. The result will be a single cell, member of the gdynprimitive class.
>> ocEP = calcep(m, X(:), 4);

% Store the equilibrium result in the model.
>> store(m, ocEP);

% Calculate the truncated time horizon based on the eigenvalues of the equilibrium.
>> epidx = 1; eigval = real(eig(ocEP{epidx})); eigval(eigval < -1e-5) = []; T = 10 / min(abs(eigval));

% Initialize the continuation process with the equilibrium point and set the continuation parameters.
>> sol = initocmat_AE_AE(m, ocEP{epidx}, 'initialstate', 'continuationcoordinate', 1:2, ...
    'targetvalue', [2.5e4 2.5e4], 'truncationtime', T);

% Start the continuation process using the 'extremal2ep' method.
>> c = bvpcont('extremal2ep', sol, [], opt);

% Store the results of the continuation process in the model.
>> store(m, 'extremal2ep');

% Retrieve the last solution path from the continuation process.
>> ocEx = extremalsolution(m);

% Save the model and results.
>> save(m, [], 1, 'eta,xi,s,h,rho,B,alpha,beta');

% Plot the solutions path and the equilibrium in the state space.
% This visualizes the extremal solution and equilibrium points.
>> clf; xcoord = 1; ycoord = 2; xvar = 'state'; yvar = 'state'; h = plotcont(m, xvar, xcoord, yvar, ycoord, ...
    'contfield', 'ExtremalSolution', 'Index', []); hold on; plotlimitset(m, xvar, xcoord, yvar, ycoord, ...
    'limitclass', 'Equilibrium', 'Index', [], 'Marker', '.', 'MarkerSize', 14, 'showspp', 1, 'showtag', 1); figure(gcf);

% Plot the time paths of the solutions and the control variables.
>> clf; xcoord = 1; ycoord = 1:3; xvar = 'time'; yvar = 'control'; h = plotcont(m, xvar, xcoord, yvar, ycoord, ...
    'contfield', 'ExtremalSolution', 'Index', []); figure(gcf);
