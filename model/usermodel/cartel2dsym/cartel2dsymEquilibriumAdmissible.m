function out=cartel2dsymEquilibriumAdmissible(t,depvar,par,arcid)
%
% B=CARTEL2DSYMEQUILIBRIUMADMISSIBLE returns a value determine the admissibilty of an equilibrium during continuation with MATCONT.
%
% this file was automatically created: 24-Aug-2024 15:01:54
% written by Dieter Grass, 2001 - 2024
	
ctrl=cartel2dsymOptimalControl(t,depvar,par,arcid);
cstr=cartel2dsymConstraint(t,depvar,par,arcid);
lagmcc=cartel2dsymLagrangeMultiplier(t,depvar,par,arcid);
	
maxfac=max(abs([ctrl;lagmcc]))^4;
switch arcid
	case 0
		out=prod(cstr([1 2 3 4]))/maxfac;
	
	case 1
		out=prod([cstr([2 3 4]);lagmcc(1)])/maxfac;
	
	case 2
		out=prod([cstr([1 3 4]);lagmcc(2)])/maxfac;
	
	case 3
		out=prod([cstr([1 2 4]);lagmcc(3)])/maxfac;
	out=cstr(4);
	case 4
		out=prod([cstr([1 2 3]);lagmcc(4)])/maxfac;
	
	case 5
		out=prod([cstr([3 4]);lagmcc([1 2])])/maxfac;
	
	case 6
		out=prod([cstr([2 4]);lagmcc([1 3])])/maxfac;
	%out=lagmcc(1);
    out=cstr(4);
	case 7
		out=prod([cstr([2 3]);lagmcc([1 4])])/maxfac;
	
	case 8
		out=prod([cstr([1 4]);lagmcc([2 3])])/maxfac;
	
	case 9
		out=prod([cstr([1 3]);lagmcc([2 4])])/maxfac;
	
	case 10
		out=prod([cstr([1 2]);lagmcc([3 4])])/maxfac;
	
	case 11
		out=prod([cstr(4);lagmcc([1 2 3])])/maxfac;
	
	case 12
		out=prod([cstr(3);lagmcc([1 2 4])])/maxfac;
	
	case 13
		out=prod([cstr(2);lagmcc([1 3 4])])/maxfac;
	
	case 14
		out=prod([cstr(1);lagmcc([2 3 4])])/maxfac;
	
end
if out<0
	out
end
