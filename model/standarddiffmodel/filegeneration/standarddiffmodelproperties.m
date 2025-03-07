function [mainsectionnames mandatorysectionnames]=standarddiffmodelproperties(mainsection)
%
%
mandatorysectionnames='';
if nargin==0
    mainsection='init';
end

switch mainsection
    case 'init'
        mainsectionnames={'variable','statedynamics','objective','arcdefinition', ...
            'control','salvagevalue','optimizationtype','controlconstraint','stateconstraint', ...
            'parameter','maximizingderivativevariable','maximizingexplicitvariable','description', ...
            'exogenousfunction','maptype'};
        mandatorysectionnames=mainsectionnames(1:3);
    case 'foc'
        mainsectionnames={'abbreviation','pontryaginfunction','pontryaginfunctionderivativeorder1','adjointsystem' ...
            'maximizingcondition','canonicalsystem','pontryaginfunctionderivativeorder2','canonicalsystemderivativeorder1'};
        mandatorysectionnames=mainsectionnames(2:5);
end
