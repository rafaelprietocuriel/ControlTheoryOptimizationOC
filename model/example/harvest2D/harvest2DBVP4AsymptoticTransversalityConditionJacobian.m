function [Ja Jb Jpar]=harvest2DBVP4AsymptoticTransversalityConditionJacobian(dynVar,LimitsetType,AsymptoticBCMatrix,SaddlePoint)
%
% the BC file for the asymptotic transversality condition
                                                                                                          
out=[];
switch LimitsetType
	case 'EP'
        Ja=[0 0 0 0;0 0 0 0];
		Jb=AsymptoticBCMatrix;
        Jpar=[0;0];
	case 'LC'
end
