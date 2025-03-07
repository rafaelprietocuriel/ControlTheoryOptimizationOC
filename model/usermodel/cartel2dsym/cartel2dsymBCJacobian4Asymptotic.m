function [Ja Jb Jpar]=cartel2dsymBCJacobian4Asymptotic(depvarb,arcid,asymptoticmatrix,saddlepoint)
%
% this file was automatically created: 24-Aug-2024 15:01:54
% written by Dieter Grass, 2001 - 2024
Ja=[];
Jb=[];
Jpar=[];
	
switch arcid
	case 0
		if size(asymptoticmatrix,1)==2
			Ja=[0,0,0,0,0; ...
				0,0,0,0,0];
			Jb=[0,0,0,0,0; ...
				0,0,0,0,0];
			Jb(:,1:4)=asymptoticmatrix(:,1:4);
		else
		end
	
	case 1
		if size(asymptoticmatrix,1)==2
			Ja=[0,0,0,0,0; ...
				0,0,0,0,0];
			Jb=[0,0,0,0,0; ...
				0,0,0,0,0];
			Jb(:,1:4)=asymptoticmatrix(:,1:4);
		else
		end
	
	case 2
		if size(asymptoticmatrix,1)==2
			Ja=[0,0,0,0,0; ...
				0,0,0,0,0];
			Jb=[0,0,0,0,0; ...
				0,0,0,0,0];
			Jb(:,1:4)=asymptoticmatrix(:,1:4);
		else
		end
	
	case 3
		if size(asymptoticmatrix,1)==2
			Ja=[0,0,0,0; ...
				0,0,0,0];
			Jb=[0,0,0,0; ...
				0,0,0,0];
			Jb(:,1:4)=asymptoticmatrix(:,1:4);
		else
		end
	
	case 4
		if size(asymptoticmatrix,1)==2
			Ja=[0,0,0,0,0,0; ...
				0,0,0,0,0,0];
			Jb=[0,0,0,0,0,0; ...
				0,0,0,0,0,0];
			Jb(:,1:4)=asymptoticmatrix(:,1:4);
		else
		end
	
	case 5
		if size(asymptoticmatrix,1)==2
			Ja=[0,0,0,0,0; ...
				0,0,0,0,0];
			Jb=[0,0,0,0,0; ...
				0,0,0,0,0];
			Jb(:,1:4)=asymptoticmatrix(:,1:4);
		else
		end
	
	case 6
		if size(asymptoticmatrix,1)==2
			Ja=[0,0,0,0; ...
				0,0,0,0];
			Jb=[0,0,0,0; ...
				0,0,0,0];
			Jb(:,1:4)=asymptoticmatrix(:,1:4);
		else
		end
	
	case 7
		if size(asymptoticmatrix,1)==2
			Ja=[0,0,0,0,0; ...
				0,0,0,0,0];
			Jb=[0,0,0,0,0; ...
				0,0,0,0,0];
			Jb(:,1:4)=asymptoticmatrix(:,1:4);
		else
		end
	
	case 8
		if size(asymptoticmatrix,1)==2
			Ja=[0,0,0,0; ...
				0,0,0,0];
			Jb=[0,0,0,0; ...
				0,0,0,0];
			Jb(:,1:4)=asymptoticmatrix(:,1:4);
		else
		end
	
	case 9
		if size(asymptoticmatrix,1)==2
			Ja=[0,0,0,0,0; ...
				0,0,0,0,0];
			Jb=[0,0,0,0,0; ...
				0,0,0,0,0];
			Jb(:,1:4)=asymptoticmatrix(:,1:4);
		else
		end
	
	case 10
		if size(asymptoticmatrix,1)==2
			Ja=[0,0,0,0,0; ...
				0,0,0,0,0];
			Jb=[0,0,0,0,0; ...
				0,0,0,0,0];
			Jb(:,1:4)=asymptoticmatrix(:,1:4);
		else
		end
	
	case 11
		if size(asymptoticmatrix,1)==2
			Ja=[0,0,0,0; ...
				0,0,0,0];
			Jb=[0,0,0,0; ...
				0,0,0,0];
			Jb(:,1:4)=asymptoticmatrix(:,1:4);
		else
		end
	
	case 12
		if size(asymptoticmatrix,1)==2
			Ja=[0,0,0,0; ...
				0,0,0,0];
			Jb=[0,0,0,0; ...
				0,0,0,0];
			Jb(:,1:4)=asymptoticmatrix(:,1:4);
		else
		end
	
	case 13
		if size(asymptoticmatrix,1)==2
			Ja=[0,0,0,0; ...
				0,0,0,0];
			Jb=[0,0,0,0; ...
				0,0,0,0];
			Jb(:,1:4)=asymptoticmatrix(:,1:4);
		else
		end
	
	case 14
		if size(asymptoticmatrix,1)==2
			Ja=[0,0,0,0; ...
				0,0,0,0];
			Jb=[0,0,0,0; ...
				0,0,0,0];
			Jb(:,1:4)=asymptoticmatrix(:,1:4);
		else
		end
	
end
