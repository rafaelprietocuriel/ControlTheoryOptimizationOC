function [mainnames,mandatorynames]=ppdemodelvariables()
%
%
mainnames={'state','control','independent','time','space','femdata','costate','spacegeometry','boundarycondition','endtime','lagrangemultcc','lagrangemultsc'};
mandatorynames=mainnames(1:4);