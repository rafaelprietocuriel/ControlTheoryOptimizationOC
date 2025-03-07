function [mainnames,mandatorynames]=standardmodelvariables()
%
%
mainnames={'state','control','independent','costate','endtime','lagrangemultcc','lagrangemultsc','connectiontime','variationparameter','variationstate','variationcostate','variationcontrol','exogenousstate'};
mandatorynames=mainnames(1:2);