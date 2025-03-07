function [Ja Jb Jpar]=harvest2DBVP4InitialStatesJacobian(dynVara,StartValue,ContinuationVector,TargetCoordinate,contparameter)
 
Ja=[0 0 0 0;0 0 0 0];
Ja(TargetCoordinate(1),TargetCoordinate(1))=1;
Ja(TargetCoordinate(2),TargetCoordinate(2))=1;
Jb=[0 0 0 0;0 0 0 0];
Jpar=-ContinuationVector;
