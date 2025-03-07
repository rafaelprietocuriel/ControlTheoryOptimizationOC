function [mainnames mandatorynames]=standarddiffmodelvariables()
%
%
mainnames={'state','control','independent','costate','endtime','lagrangemultcc','lagrangemultsc'};
mandatorynames=mainnames(1:2);