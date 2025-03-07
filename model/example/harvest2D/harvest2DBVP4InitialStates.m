function out = harvest2DBVP4InitialStates(dynVara,StartValue,ContinuationVector,TargetCoordinate,contparameter)
%
% the BC file returning the continuation residuum
                                                                                                                                                                                                      
% global variable
                                                                                                                                                                                                      
out=dynVara(TargetCoordinate)-StartValue-contparameter*ContinuationVector;
