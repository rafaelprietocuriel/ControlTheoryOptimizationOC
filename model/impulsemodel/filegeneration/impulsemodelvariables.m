function [mainnames mandatorynames]=impulsemodelvariables()
%
%
mainnames={'state','control','icontrol','independent','costate','impulsetime','endtime','lagrangemultcc','lagrangemultsc'};
mandatorynames=mainnames(1:2);