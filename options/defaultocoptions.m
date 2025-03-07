function opt = defaultocoptions()
%
% DEFAULTOCOPTIONS Creates and returns a structure containing all default
% options for the OCMAT toolbox.
%
% The intention of this structure is to allow a compact handling of the
% different options. It consists of several substructures each containing
% default options needed for different computations of the toolbox.
%
% OPT = DEFAULTOCOPTIONS the structure OPT consists of several
% substructures, which can be identified by their name and which contain
% options needed for different applications within the toolbox: 
%       'OC'        .....   general OCMat options
%       'OCCONT'    .....   options for the continuation algorithm
%       'OCCONTARG' .....   options needed for the initialization of a BVP
%       'ODE'       .....   options for solving an IVP (provided by MATLAB
%               - see odeset / odeget) 
%       'BVP'       .....   options for solving an IVP (provided by MATLAB
%               - see bvpset / bvpget) 
%       'EP'        .....   options needed for the calculation of steady
%                           states (provided by MATLAB (optimization
%                           toolbox)- see optimset/optimget)
%       'MATCONT'   .....   contains options required by MATCONT toolbox
%               (see contset / contget) 
%       'SBVP'      .....   for the BVPSuite (not implemented yet)
%       'GRADIENT'  .....   options for the gradient algorithm
%
% Note that the default values of the options are certainly not perfectly
% suitable for all models. However, the options can be easily adjusted by
% using the function 'setocoptions', by which one can access and manipulate
% the values of this option structure. Please be aware that a badly done
% (manual) manipulation or the use of a different option structure
% might lead to a malfunction of the toolbox.

%matlabversion=ver;

optimizationToolboxFlag=1;%~isempty(strfind([matlabversion.Name],'Optimization Toolbox'));
% OC Continuation Input Argument Options - Needed for the
% initialization of a structure specifying the BVP and some default values for the
% continuation algorithm
defaultoccontarg=struct('Increment', 1e-5,'MaxCorrIters',10,'MaxTestIters',10,'FunTolerance',1e-6,'VarTolerance',1e-6,'MaxContinuationSteps',300, ...
    'TestTolerance',1e-5,'Singularities',0,'Backward',0,'InitStepWidth',0.01,'MaxStepWidth',0.1, ...
    'MinStepWidth', 1e-5,'IncreaseFactor',1.3,'DecreaseFactor',0.5,'ContinuationMethod',1,'WorkSpace',0, ...
    'HitTargetValue',1,'CheckAngle',0.9,'CheckStep',3,'SaveIntermediate','on','PrintContStats','on','PlotCont','on', ...
    'IgnoreSingularity',[],'ExitOnTargetValue',1,'Locators',[],'Userfunctions',0,'CheckAdmissibility',1,'Adapt',3, ...
    'TotalRelativeDistance',1,'MaxDistance',inf,'MaxDistanceCoordinate',[],'AsymptoticBCMethod','schur');

% NEWTON Input Argument Options - 
defaultocnewton=struct('MaxNewtonIters',4,'MaxProbes',4,'AbsTol',1e-6,'RelTol',1e-3,'Display','', ...
    'TRM',0,'Log',0,'LambdaMin',0.001,'UpdateJacFactor',0.5,'SwitchToFFNFactor',0.5,'CheckSingular',1,'SingularThreshold',1e-1, ...
    'GeneralizedDerivativeThreshold',1e-6);

defaultocgeneral=struct('BVPMethod','bvp4c','ODESolver','ode45','EquationSolver','fsolve', ...
    'GradMethod','fmincon','ZeroTimeDiffTolerance',1e-6,'ZeroDeviationTolerance',1e-6,'AdmissibleTolerance',1e-6,'ImaginaryTolerance',1e-10, ...
    'TrivialArcMeshNum',5,'ContLog',0,'NewtonSolver','newtcorr4bvp','LegendreClebsch','off','ShootMethod','ocshootloc');

% ODE Options - Used for the solution of IVP problems - See MATLAB help
defaultode=odeset('AbsTol',1e-5,'BDF','off','JConstant','off',...
    'MassSingular','maybe','NormControl','off','RelTol',1e-8,...
    'Stats','off','Vectorized','on');

%Initialization Options - Used for initializing an ocmat problem
defaultinit=struct('Simplify','off','Jacobian','explicit','Hessian','explicit','ControlDynamics','explicit', ...
    'ParameterJacobian','explicit','Force','off','MessageDisplay','off', ...
    'TestConsistency','on','TestModus','off');

defaultoptimset=optimset('GradObj','on');
%BVP Options - Used for the solution of BVP problems - See MATLAB help
try
    defaultbvp=bvpset('AbsTol',1e-6,'RelTol',1e-3,'Stats','off','Vectorized','on','NMax',4000);
    % MaxNewPts is only a field in the adapted bvpset for bvp6c
catch
    defaultbvp=bvpset('AbsTol',1e-6,'RelTol',1e-3,'Stats','off','Vectorized','on','NMax',4000);
end
% add fields for BVP Solver TOM
defaultbvp.Itnlmax=6;
defaultbvp.Monitor=1;
defaultbvp.Order=6;
defaultbvp.LUsw=1;
defaultbvp.MassMatrix=[];


%GRADIENT Options
% defaultgradient =struct('ConditionNumber',1,'LipschitzConst',1e4,'MaxIter',1e3,'Display','iter','LineSrch','ValCubic', ...
%     'Projection',0,'ODESolver','heunloc','PenaltyParameter',10,'InexactStepSizeParameter',0.5, ...
%     'GradientTolerance',1e-5,'SearchDirectionParameter',0.5,'MaxGradientIteration',1000,'MaxLinesearchIteration',30,'InitLineStepWidth',1e-2,'MinLineStepWidth',1e-10,'MaxLineStepWidth',1e2,'QPSubProblemNumeric',0, ...
%     'DirectionalDerivativeStep',1e-5,'MakeMovie',0,'GradientMappingMethod','kkt','GradientMappingGamma',30,'DirectionCorrectionType','', ...
%     'LineStepIncreaseFactor',1.5,'ObjectiveJumpThreshold',1e2,'NewtonSolver','newtcorr4grad');
defaultgradient =struct('Projection',0,'PenaltyParameter',10,'InexactStepSizeParameter',0.5,'GradientTolerance',1e-4, ...
    'SearchDirectionParameter',0.5,'MaxGradientIteration',1000,'MaxLinesearchIteration',30,'InitLineStepWidth',1e-2, ...
    'MinLineStepWidth',1e-10,'MaxLineStepWidth',1e2,'DirectionalDerivativeStep',1e-5,'MakeMovie',0,'GradientMappingMethod','nesterov', ...
    'GradientMappingGamma',30,'LineStepIncreaseFactor',1.5,'GlobalCorrector',0,'PrintGradientInfo',0,'LineSearchMethod','min');

% Equilibria Options - Used e.g. for the calculation of Equilibria - See
% MATLAB help
if optimizationToolboxFlag
    defaulteq=optimset(optimset,'Display','off','Jacobian','off','MaxFunEvals',5000,'MaxIter',5000,'TolFun',1e-10);
else
    defaulteq=struct('Display','off','Jacobian','off');
end

% Matcont Options - Used by MATCONT - See MATCONT manual
if exist('contset','file')
    defaultmatcont=contset;
    defaultmatcont=contset(defaultmatcont,'Singularities',1);
    defaultmatcont.PlotCont=1;
    defaultmatcont.SaveIntermediate='on';
else
    defaultmatcont=struct;
end

%SBVPOC Options
defaultsbvpoc =struct('MeshAdaptation',10,'MeshAdaptAbsTol',1e-6,'MeshAdaptRelTol',1e-3,'MeshAdaptFactorMax',1.5, ...
    'MeshAdaptFineMesh',1,'Vectorized','on','NMax',1000,'FJacobian',1,'BCJacobian',0,'ICJacobian',1, ...
    'CollocationNumber',3,'CollocationMethod','uniform','InterpolationMethod','spline','MeshAdaptMaxIter',18, ...
    'MeshUpdateMode',1,'KRange',100);

if optimizationToolboxFlag
    defaultzfopt=optimset('Display','off','MaxIter',90000000,'MaxFunEvals',90000000);
else
    defaultzfopt=[];
end
%Option structure
opt=struct('OCCONTARG',defaultoccontarg,'ODE',defaultode,'BVP',defaultbvp,'EQ',defaulteq, ...
    'MATCONT',defaultmatcont,'SBVPOC',defaultsbvpoc,'INIT',defaultinit,'NEWTON',defaultocnewton, ...
    'GENERAL',defaultocgeneral,'GRADIENT',defaultgradient,'STATICOPTIM',defaultoptimset);



