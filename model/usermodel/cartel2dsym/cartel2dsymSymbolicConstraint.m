function out=cartel2dsymSymbolicConstraint(arcid)
% returns the symbolic value of the exogenously given functions
%
% this file was automatically created: 24-Aug-2024 15:01:54
% written by Dieter Grass, 2001 - 2024
	
out=mystr2sym(['[u1-(ulow);' ...
	'u2-(ulow);' ...
	'v-(vlow);' ...
	'B-(u1+u2+v)]']);
