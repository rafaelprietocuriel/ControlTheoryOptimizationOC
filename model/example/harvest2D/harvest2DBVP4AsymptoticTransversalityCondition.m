function out = harvest2DBVP4AsymptoticTransversalityCondition(dynVar,LimitsetType,AsymptoticBCMatrix,SaddlePoint)
%
% the BC file for the asymptotic transversality condition
                                                                                                          
out=[];
switch LimitsetType
	case 'EP'
		out=AsymptoticBCMatrix*(dynVar(1:4)-SaddlePoint);
		%if ~isempty(OcBVPVar.ConstantCoordinates)
			%out=[out;dynVar(OcBVPVar.ConstantCoordinates)-SaddlePoint(OcBVPVar.ConstantCoordinates)];
		%end
	case 'LC'
end
