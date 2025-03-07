function [Ja Jb Jpar]=cartel2dsymBCJacobian4Initial(depvara,freepar,arcid,targetcoordinate,continuationvector)
% returns the Jacobian of the initial boundary conditions for the saddle path continuation
%
% this file was automatically created: 24-Aug-2024 15:01:54
% written by Dieter Grass, 2001 - 2024
	
Ja=[];
Jb=[];
Jpar=[];
	
switch arcid
	case 0
		if numel(targetcoordinate)==2
			Ja=[0,0,0,0,0; ...
				0,0,0,0,0];
			Jb=[0,0,0,0,0; ...
				0,0,0,0,0];
			Ja(targetcoordinate(1),targetcoordinate(1))=1;
			Ja(targetcoordinate(2),targetcoordinate(2))=1;
		else
		end
	
	case 1
		if numel(targetcoordinate)==2
			Ja=[0,0,0,0,0; ...
				0,0,0,0,0];
			Jb=[0,0,0,0,0; ...
				0,0,0,0,0];
			Ja(targetcoordinate(1),targetcoordinate(1))=1;
			Ja(targetcoordinate(2),targetcoordinate(2))=1;
		else
		end
	
	case 2
		if numel(targetcoordinate)==2
			Ja=[0,0,0,0,0; ...
				0,0,0,0,0];
			Jb=[0,0,0,0,0; ...
				0,0,0,0,0];
			Ja(targetcoordinate(1),targetcoordinate(1))=1;
			Ja(targetcoordinate(2),targetcoordinate(2))=1;
		else
		end
	
	case 3
		if numel(targetcoordinate)==2
			Ja=[0,0,0,0; ...
				0,0,0,0];
			Jb=[0,0,0,0; ...
				0,0,0,0];
			Ja(targetcoordinate(1),targetcoordinate(1))=1;
			Ja(targetcoordinate(2),targetcoordinate(2))=1;
		else
		end
	
	case 4
		if numel(targetcoordinate)==2
			Ja=[0,0,0,0,0,0; ...
				0,0,0,0,0,0];
			Jb=[0,0,0,0,0,0; ...
				0,0,0,0,0,0];
			Ja(targetcoordinate(1),targetcoordinate(1))=1;
			Ja(targetcoordinate(2),targetcoordinate(2))=1;
		else
		end
	
	case 5
		if numel(targetcoordinate)==2
			Ja=[0,0,0,0,0; ...
				0,0,0,0,0];
			Jb=[0,0,0,0,0; ...
				0,0,0,0,0];
			Ja(targetcoordinate(1),targetcoordinate(1))=1;
			Ja(targetcoordinate(2),targetcoordinate(2))=1;
		else
		end
	
	case 6
		if numel(targetcoordinate)==2
			Ja=[0,0,0,0; ...
				0,0,0,0];
			Jb=[0,0,0,0; ...
				0,0,0,0];
			Ja(targetcoordinate(1),targetcoordinate(1))=1;
			Ja(targetcoordinate(2),targetcoordinate(2))=1;
		else
		end
	
	case 7
		if numel(targetcoordinate)==2
			Ja=[0,0,0,0,0; ...
				0,0,0,0,0];
			Jb=[0,0,0,0,0; ...
				0,0,0,0,0];
			Ja(targetcoordinate(1),targetcoordinate(1))=1;
			Ja(targetcoordinate(2),targetcoordinate(2))=1;
		else
		end
	
	case 8
		if numel(targetcoordinate)==2
			Ja=[0,0,0,0; ...
				0,0,0,0];
			Jb=[0,0,0,0; ...
				0,0,0,0];
			Ja(targetcoordinate(1),targetcoordinate(1))=1;
			Ja(targetcoordinate(2),targetcoordinate(2))=1;
		else
		end
	
	case 9
		if numel(targetcoordinate)==2
			Ja=[0,0,0,0,0; ...
				0,0,0,0,0];
			Jb=[0,0,0,0,0; ...
				0,0,0,0,0];
			Ja(targetcoordinate(1),targetcoordinate(1))=1;
			Ja(targetcoordinate(2),targetcoordinate(2))=1;
		else
		end
	
	case 10
		if numel(targetcoordinate)==2
			Ja=[0,0,0,0,0; ...
				0,0,0,0,0];
			Jb=[0,0,0,0,0; ...
				0,0,0,0,0];
			Ja(targetcoordinate(1),targetcoordinate(1))=1;
			Ja(targetcoordinate(2),targetcoordinate(2))=1;
		else
		end
	
	case 11
		if numel(targetcoordinate)==2
			Ja=[0,0,0,0; ...
				0,0,0,0];
			Jb=[0,0,0,0; ...
				0,0,0,0];
			Ja(targetcoordinate(1),targetcoordinate(1))=1;
			Ja(targetcoordinate(2),targetcoordinate(2))=1;
		else
		end
	
	case 12
		if numel(targetcoordinate)==2
			Ja=[0,0,0,0; ...
				0,0,0,0];
			Jb=[0,0,0,0; ...
				0,0,0,0];
			Ja(targetcoordinate(1),targetcoordinate(1))=1;
			Ja(targetcoordinate(2),targetcoordinate(2))=1;
		else
		end
	
	case 13
		if numel(targetcoordinate)==2
			Ja=[0,0,0,0; ...
				0,0,0,0];
			Jb=[0,0,0,0; ...
				0,0,0,0];
			Ja(targetcoordinate(1),targetcoordinate(1))=1;
			Ja(targetcoordinate(2),targetcoordinate(2))=1;
		else
		end
	
	case 14
		if numel(targetcoordinate)==2
			Ja=[0,0,0,0; ...
				0,0,0,0];
			Jb=[0,0,0,0; ...
				0,0,0,0];
			Ja(targetcoordinate(1),targetcoordinate(1))=1;
			Ja(targetcoordinate(2),targetcoordinate(2))=1;
		else
		end
	
end
if numel(targetcoordinate)==2
	numfree=numel(freepar);
	Jpar=zeros(2,numfree);
	Jpar(targetcoordinate,numfree)=-continuationvector;
end
