function h=cartel2dsymPlotIndifferenceContinuation(t,depvar,par,arcid,freepar,tangent,plotflag)
% plotting commands called during the continuation process
% the user can adapt this file to her/his needs
%
% this file was automatically created: 24-Aug-2024 15:01:54
% written by Dieter Grass, 2001 - 2024
	
global OCMATCONT OCMATINDIF
	
color=lines;
numarc=numel(arcid);
leftarcindex=OCMATCONT.HE.TIMEDDATA.leftarcindex;
rightarcindex=OCMATCONT.HE.TIMEDDATA.rightarcindex;
h=zeros(1,numarc);
	
for ii=1:numarc
	h(ii)=plot(depvar(1,leftarcindex(ii):rightarcindex(ii)),depvar(2,leftarcindex(ii):rightarcindex(ii)));
	set(h(ii),'Color',color(arcid(ii)+1,:));
	if ii==1
		hold on
	end
end
%set(h(OCMATINDIF.arccoord{OCMATINDIF.indifferenceorder}),'LineStyle','--')
hold off
	
xlabel('$C_1$','Interpreter','Latex')
ylabel('$C_2$','Interpreter','Latex')
