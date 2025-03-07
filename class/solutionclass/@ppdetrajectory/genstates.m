function out=genstates(ppdeTrj)

genstatesidx=ppdeTrj.discretizationinfo.coord.totaldepvar;
out=ppdeTrj.y(genstatesidx,:);