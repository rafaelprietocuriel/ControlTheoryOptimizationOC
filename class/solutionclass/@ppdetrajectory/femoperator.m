function varargout=femoperator(ppdeTrj)

varargout{1}=ppdeTrj.discretizationinfo.femop.M;
varargout{2}=ppdeTrj.discretizationinfo.femop.invM;
varargout{3}=ppdeTrj.discretizationinfo.femop.bcG;
varargout{4}=ppdeTrj.discretizationinfo.femop.K;
varargout{5}=ppdeTrj.discretizationinfo.femop.Kadv;